#include "vcf-io.h"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <filesystem>  // NOLINT(build/c++17)
#include <memory>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <xoos/io/alignment-reader.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/htslib-util/read-util.h>
#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "compute-vcf-features/vcf-header-util.h"
#include "phaser/params.h"
#include "phaser/variant.h"

namespace xoos::svc::phaser {

// A phase block needs at least two variants before it can be merged into an MNV.
static constexpr size_t kMinBlockSize = 2;

VcfData ReadVcf(const fs::path& vcf_path, const Params& params) {
  VcfData data;

  const io::VcfReader reader(vcf_path);
  data.header = reader.GetHeader();

  const s32 n_samples = data.header->GetNumSamples();
  // Resolve the tumor sample column from the VCF header metadata (`tumor_sample`), so we annotate
  // the correct sample regardless of column order or a tumor-only VCF. When the header lacks the
  // metadata, fall back to the last sample (tumor is the trailing column in the standard T/N layout,
  // and the sole column in a tumor-only VCF). Guard against out-of-range indices so downstream
  // FORMAT array indexing stays in bounds.
  s32 resolved_col = n_samples - 1;
  if (const auto tn_idx = GetTumorNormalSampleIndexes(data.header); tn_idx.has_value()) {
    resolved_col = tn_idx->tumor_sample_idx;
  }
  data.tumor_sample_idx = (resolved_col >= 0 && resolved_col < n_samples) ? resolved_col : 0;

  data.records = reader.GetAllRecords();

  // A bi-allelic site has exactly one REF and one ALT allele.
  static constexpr s32 kBiallelicAlleleCount = 2;

  s32 filtered_count = 0;
  // SonarQube misidentifies `data` as the loop-control variable because `data.records.size()` appears in the
  // loop condition; the actual counter `rec_idx` is never modified in the body, so the rule does not apply.
  for (s32 rec_idx = 0; rec_idx < static_cast<s32>(data.records.size()); ++rec_idx) {  // NOSONAR(cpp:S892)
    const auto& record = data.records[static_cast<size_t>(rec_idx)];
    // Only bi-allelic sites are considered (one REF, one ALT).
    if (record->NumAlleles() != kBiallelicAlleleCount) {
      continue;
    }

    const bool is_pass = record->IsPass();
    if (!is_pass && !params.include_filtered && !params.rescue_filtered) {
      ++filtered_count;
      continue;
    }

    Variant v;
    v.chrom = record->Chromosome();
    // 0-based
    v.pos = record->Position();
    v.ref = record->Allele(0);
    v.alt = record->Allele(1);
    v.vcf_record_index = rec_idx;
    v.is_pass = is_pass;

    const auto filters = record->GetFilters();
    v.filter = filters.empty() ? "." : fmt::format("{}", fmt::join(filters, ";"));

    data.variants.push_back(std::move(v));
  }

  Logging::Info("Read {} variants from VCF ({} filtered out)", data.variants.size(), filtered_count);
  return data;
}

namespace {

// Determine if a phase block contains only SNVs.
bool IsSnvOnlyBlock(const vec<s32>& sorted_indices, const vec<Variant>& variants) {
  return std::ranges::all_of(sorted_indices,
                             [&variants](const s32 vi) { return variants[static_cast<size_t>(vi)].IsSnv(); });
}

// Extract reference bases from a BAM read at given reference positions.
// Walks the CIGAR to map reference positions to query bases.
// Returns the reference-aligned bases (matches only, skipping insertions).
std::string ExtractRefBasesFromRead(bam1_t* const aln, const s64 ref_start, const s64 ref_end) {
  const u8* const seq = bam_get_seq(aln);
  const std::span<const u32> cigar{bam_get_cigar(aln), aln->core.n_cigar};

  std::string result(static_cast<size_t>(ref_end - ref_start), 'N');
  s64 rpos = aln->core.pos;
  s64 qpos = 0;

  for (const u32 cigar_op : cigar) {
    const auto op = static_cast<s32>(bam_cigar_op(cigar_op));
    const auto len = static_cast<s32>(bam_cigar_oplen(cigar_op));

    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
      for (s32 j = 0; j < len; ++j) {
        if (rpos >= ref_start && rpos < ref_end) {
          result[static_cast<size_t>(rpos - ref_start)] = io::GetBase(seq, static_cast<u32>(qpos));
        }
        ++rpos;
        ++qpos;
      }
    } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
      qpos += len;
    } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
      rpos += len;
    }
    if (rpos >= ref_end) {
      break;
    }
  }
  return result;
}

// Fetch reference sequence for a region from the BAM by finding a spanning read.
// Falls back to 'N' if no suitable read is found.
std::string FetchRefFromBam(const io::AlignmentReader& bam_handle,
                            const std::string& chrom,
                            const s64 ref_start,
                            const s64 ref_end) {
  // A run of N of the requested length, used when the reference cannot be resolved.
  const auto n_fill = [ref_start, ref_end] { return std::string(static_cast<size_t>(ref_end - ref_start), 'N'); };

  if (bam_handle.bam_file == nullptr || bam_handle.hdr == nullptr || bam_handle.idx == nullptr) {
    return n_fill();
  }

  const s32 tid = sam_hdr_name2tid(bam_handle.hdr.get(), chrom.c_str());
  if (tid < 0) {
    return n_fill();
  }

  // unique_ptr with custom deleter adopts the htslib C-API pointer; make_unique cannot be used.
  const io::HtsItrPtr iter{sam_itr_queryi(bam_handle.idx.get(), tid, ref_start, ref_end)};  // NOSONAR(cpp:M23_286)
  std::string result;

  if (iter != nullptr) {
    const io::Bam1Ptr aln_ptr = io::BamInit1();
    bam1_t* const aln = aln_ptr.get();
    while (sam_itr_next(bam_handle.bam_file.get(), iter.get(), aln) >= 0) {
      // htslib BAM_F* flag constants are signed int macros; cast to unsigned so the mask is a
      // clean unsigned bitwise op (matches compute-bam-region-features.cpp).
      static constexpr auto kSkipMask =
          static_cast<u32>(BAM_FUNMAP) | static_cast<u32>(BAM_FSECONDARY) | static_cast<u32>(BAM_FSUPPLEMENTARY);
      if ((static_cast<u32>(aln->core.flag) & kSkipMask) != 0) {
        continue;
      }
      // Need a read that spans the entire region.
      if (aln->core.pos <= ref_start && bam_endpos(aln) >= ref_end) {
        result = ExtractRefBasesFromRead(aln, ref_start, ref_end);
        if (!result.contains('N')) {
          // Found a clean spanning read.
          break;
        }
      }
    }
  }

  if (result.empty()) {
    return n_fill();
  }
  return result;
}

// Build the reference backbone spanning [block_start, block_end) for a phase block.
// When an indexed FASTA is supplied, it is authoritative because a phase block can span positions
// that no single read covers. The BAM path remains a fallback for callers without a reference.
// Returns an empty string when gaps cannot be resolved.
std::string BuildRefBackbone(const vec<s32>& sorted_indices,
                             const vec<Variant>& variants,
                             io::FastaReader* const fasta_reader,
                             const io::AlignmentReader& bam_handle,
                             const std::string& chrom,
                             const s64 block_start,
                             const s64 block_end) {
  if (fasta_reader != nullptr) {
    try {
      const std::string reference = fasta_reader->GetSequence(chrom, block_start, block_end);
      if (reference.size() == static_cast<size_t>(block_end - block_start)) {
        return reference;
      }
      Logging::Warn("Reference sequence length mismatch for {}:{}-{}", chrom, block_start + 1, block_end);
    } catch (const std::exception& e) {  // NOSONAR(cpp:S1181) — reference fetch failure has a safe BAM fallback
      Logging::Warn(
          "Failed to fetch reference sequence for {}:{}-{} ({})", chrom, block_start + 1, block_end, e.what());
    }
  }

  std::string ref_backbone(static_cast<size_t>(block_end - block_start), 'N');
  for (const s32 vi : sorted_indices) {
    const auto& v = variants[static_cast<size_t>(vi)];
    const auto offset = static_cast<size_t>(v.pos - block_start);
    for (size_t i = 0; i < v.ref.size() && (offset + i) < ref_backbone.size(); ++i) {
      ref_backbone[offset + i] = v.ref[i];
    }
  }

  // Fill remaining gap positions from the BAM, if one is available.
  if (ref_backbone.contains('N') && bam_handle.bam_file != nullptr) {
    const std::string bam_ref = FetchRefFromBam(bam_handle, chrom, block_start, block_end);
    for (size_t i = 0; i < ref_backbone.size() && i < bam_ref.size(); ++i) {
      if (ref_backbone[i] == 'N') {
        ref_backbone[i] = bam_ref[i];
      }
    }
  }

  // If gaps remain, the reference could not be resolved.
  if (ref_backbone.contains('N')) {
    return {};
  }
  return ref_backbone;
}

// Append the contribution of a single overlapping variant (v.pos < cursor) to merged_alt,
// advancing cursor. Returns true when the variant contributed alt bases.
bool AppendOverlappingVariant(const Variant& v, std::string& merged_alt, s64& cursor) {
  if (v.EndPos() <= cursor) {
    // A deletion can completely cover the reference anchor of a later insertion.
    // In that case the anchor itself is already represented by the deletion, but
    // the insertion suffix still contributes to the combined alternate allele.
    // For example, TCCA>T plus A>ATT becomes TCCA>TTT.
    if (v.IsInsertion() && v.alt.starts_with(v.ref)) {
      merged_alt += v.alt.substr(v.ref.size());
      return true;
    }
    return false;
  }
  const auto skip = static_cast<size_t>(cursor - v.pos);
  if (skip >= v.alt.size() || skip >= v.ref.size()) {
    return false;
  }
  merged_alt += v.alt.substr(skip);
  cursor = v.EndPos();
  return true;
}

// Assemble the merged ALT allele by walking the sorted variants across the backbone.
// Returns an empty string when the layout cannot be cleanly merged.
std::string BuildMergedAlt(const vec<s32>& sorted_indices,
                           const vec<Variant>& variants,
                           const std::string_view ref_backbone,
                           const s64 block_start,
                           const s64 block_end) {
  std::string merged_alt;
  s64 cursor = block_start;

  for (const s32 vi : sorted_indices) {
    const auto& v = variants[static_cast<size_t>(vi)];
    if (v.pos < cursor) {
      // Overlapping variant (e.g. deletion + insertion). If its contribution cannot be cleanly
      // appended, abort merging for the whole block rather than emitting a merged ALT that drops
      // this variant's alleles — the caller (MergeBlock) skips the block on an empty result.
      if (!AppendOverlappingVariant(v, merged_alt, cursor)) {
        return {};
      }
      continue;
    }
    // Append reference between cursor and this variant.
    if (v.pos > cursor) {
      const auto ref_off = static_cast<size_t>(cursor - block_start);
      const auto ref_len = static_cast<size_t>(v.pos - cursor);
      if (ref_off + ref_len > ref_backbone.size()) {
        return {};
      }
      merged_alt += ref_backbone.substr(ref_off, ref_len);
    }
    merged_alt += v.alt;
    cursor = v.EndPos();
  }

  // Append any trailing reference.
  if (cursor < block_end) {
    const auto ref_off = static_cast<size_t>(cursor - block_start);
    if (ref_off <= ref_backbone.size()) {
      merged_alt += ref_backbone.substr(ref_off);
    }
  }
  return merged_alt;
}

// Build merged REF and ALT alleles for a phase block.
// For SNV-only blocks: fills gaps between SNVs with reference bases.
// For complex blocks (with indels): builds the combined haplotype.
// Returns {merged_ref, merged_alt} or empty strings if merging is not possible.
std::pair<std::string, std::string> BuildMergedAlleles(const vec<s32>& sorted_indices,
                                                       const vec<Variant>& variants,
                                                       io::FastaReader* const fasta_reader,
                                                       const io::AlignmentReader& bam_handle) {
  if (sorted_indices.size() < kMinBlockSize) {
    return {"", ""};
  }

  const auto& first = variants[static_cast<size_t>(sorted_indices.front())];
  const std::string& chrom = first.chrom;

  const s64 block_start = first.pos;
  s64 block_end = first.EndPos();
  for (const s32 vi : sorted_indices) {
    block_end = std::max(block_end, variants[static_cast<size_t>(vi)].EndPos());
  }

  const std::string merged_ref =
      BuildRefBackbone(sorted_indices, variants, fasta_reader, bam_handle, chrom, block_start, block_end);
  if (merged_ref.empty()) {
    return {"", ""};
  }

  const std::string merged_alt = BuildMergedAlt(sorted_indices, variants, merged_ref, block_start, block_end);
  if (merged_ref == merged_alt || merged_alt.empty()) {
    return {"", ""};
  }
  return {merged_ref, merged_alt};
}

// The per-block phasing annotation applied to every record in a phase set.
// Only the record and its haplotype vary between calls; these fields are constant across the block.
struct PhasingAnnotation {
  s32 n_samples;
  s32 tumor_idx;
  s32 phase_set;
  std::string_view pid;
  std::string_view mnv_tag;
  bool report_mnv;
};

// Apply phasing FORMAT fields (PGT/PID/PS) to a record for the tumor sample.
// haplotype: 0 = PGT 0|1 (default/hap A), 1 = PGT 1|0 (hap B).
void AnnotatePhasing(const io::VcfRecordPtr& record, const PhasingAnnotation& annotation, const s32 haplotype) {
  vec<std::string> pgt_vals(static_cast<size_t>(annotation.n_samples), ".");
  pgt_vals[static_cast<size_t>(annotation.tumor_idx)] = (haplotype == 0) ? "0|1" : "1|0";
  record->SetFormatField("PGT", pgt_vals);

  vec<std::string> pid_vals(static_cast<size_t>(annotation.n_samples), ".");
  pid_vals[static_cast<size_t>(annotation.tumor_idx)] = std::string{annotation.pid};
  record->SetFormatField("PID", pid_vals);

  vec<s32> ps_vals(static_cast<size_t>(annotation.n_samples), bcf_int32_missing);
  ps_vals[static_cast<size_t>(annotation.tumor_idx)] = annotation.phase_set;
  record->SetFormatField("PS", ps_vals);

  if (annotation.report_mnv && !annotation.mnv_tag.empty()) {
    record->SetInfoField("MNVTAG", vec<std::string>{std::string{annotation.mnv_tag}});
  }
}

// Toggles that shape the output header and per-block annotation.
// merge_enabled: MNV/complex merging is active; report_mnv: emit the MNVTAG INFO field.
struct PhasingOutputOptions {
  bool merge_enabled;
  bool report_mnv;
};

// Clone the input header and append the phasing FORMAT/INFO/FILTER lines.
io::VcfHeaderPtr BuildPhasedHeader(const VcfData& vcf_data, const Params& params, const PhasingOutputOptions& options) {
  const io::VcfHeaderPtr hdr = vcf_data.header->Clone();
  hdr->AddFormatLine({.id = "PGT",
                      .description = "Physical phasing haplotype",
                      .number = io::kNumberOne,
                      .type = io::FieldType::kString});
  hdr->AddFormatLine(
      {.id = "PID", .description = "Physical phasing ID", .number = io::kNumberOne, .type = io::FieldType::kString});
  hdr->AddFormatLine(
      {.id = "PS", .description = "Phase set identifier", .number = io::kNumberOne, .type = io::FieldType::kInteger});
  if (options.report_mnv) {
    hdr->AddInfoLine({.id = "MNVTAG",
                      .description = "MNV representation of phased variants",
                      .number = io::kNumberOne,
                      .type = io::FieldType::kString});
  }
  if (options.merge_enabled && params.keep_components) {
    hdr->AddFilterLine(
        {.id = "PHASED_COMPONENT",
         .description = "Original component of a phased variant; the merged representation is emitted separately as "
                        "PASS"});
  }
  if (params.command_line) {
    hdr->AddCustomMetaDataLine(io::GetCommandInfo(*params.command_line));
  }
  hdr->Sync();
  return hdr;
}

// Order a block's variant indices by position, keeping haplotypes parallel.
void OrderBlockByPosition(const PhaseBlock& block,
                          const vec<Variant>& variants,
                          vec<s32>& ordered_indices,
                          vec<s32>& ordered_haps) {
  const vec<s32>& src = block.variant_indices;
  vec<size_t> order(src.size());
  std::ranges::iota(order, size_t{0});
  std::ranges::sort(order, [&variants, &src](const size_t a, const size_t b) {
    return variants[static_cast<size_t>(src[a])].pos < variants[static_cast<size_t>(src[b])].pos;
  });
  ordered_indices.reserve(order.size());
  for (const size_t idx : order) {
    ordered_indices.push_back(src[idx]);
    if (!block.haplotypes.empty()) {
      ordered_haps.push_back(block.haplotypes[idx]);
    }
  }
}

// Shared state threaded through per-block processing.
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members) — short-lived context bound to caller-owned data
struct WriteContext {
  const VcfData& vcf_data;
  const vec<Variant>& variants;
  const Params& params;
  const io::VcfHeaderPtr& hdr;
  vec<io::VcfRecordPtr>& out_records;
  std::unordered_set<s32>& merged_away;
  std::unordered_map<s32, io::VcfRecordPtr>& merged_by_first_idx;
  io::FastaReader* fasta_reader;
  io::AlignmentReader& bam_handle;
  s32 n_samples;
  s32 tumor_idx;
  bool merge_enabled;
  bool report_mnv;
};

// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

// Build the MNV tag summarizing a block's constituent variants.
std::string BuildMnvTag(const vec<s32>& sorted_indices, const vec<Variant>& variants) {
  vec<std::string> parts;
  parts.reserve(sorted_indices.size());
  for (const s32 vi : sorted_indices) {
    const auto& v = variants[static_cast<size_t>(vi)];
    parts.push_back(fmt::format("{}:{}{}>{}", v.chrom, v.pos + 1, v.ref, v.alt));
  }
  return fmt::format("{}", fmt::join(parts, ","));
}

// Construct and register the merged MNV record for a block, if it can be built.
void MergeBlock(WriteContext& ctx,
                const vec<s32>& sorted_indices,
                const vec<s32>& ordered_haps,
                const Variant& first_var,
                const s32 ps,
                const std::string& pid,
                const std::string& mnv_tag) {
  const auto [merged_ref, merged_alt] =
      BuildMergedAlleles(sorted_indices, ctx.variants, ctx.fasta_reader, ctx.bam_handle);
  if (merged_ref.empty() || merged_alt.empty()) {
    Logging::Debug("Could not merge block at {}:{}", first_var.chrom, first_var.pos + 1);
    return;
  }

  // Clone the first constituent record as the template for the merged record.
  const s32 template_idx = ctx.variants[static_cast<size_t>(sorted_indices.front())].vcf_record_index;
  io::VcfRecordPtr merged = ctx.vcf_data.records[static_cast<size_t>(template_idx)]->Clone(ctx.hdr);
  merged->SetPosition(static_cast<s32>(first_var.pos));
  merged->SetAlleles({merged_ref, merged_alt});
  // Merged record is the canonical representation — always PASS.
  merged->SetFilter("PASS");

  // QUAL of the merged record is the minimum QUAL across its constituents.
  std::optional<f32> min_qual;
  for (const s32 vi : sorted_indices) {
    const s32 comp_idx = ctx.variants[static_cast<size_t>(vi)].vcf_record_index;
    const auto comp_qual = ctx.vcf_data.records[static_cast<size_t>(comp_idx)]->GetQuality();
    if (comp_qual.has_value() && (!min_qual.has_value() || *comp_qual < *min_qual)) {
      min_qual = comp_qual;
    }
  }
  merged->SetQuality(min_qual);
  const PhasingAnnotation annotation{.n_samples = ctx.n_samples,
                                     .tumor_idx = ctx.tumor_idx,
                                     .phase_set = ps,
                                     .pid = pid,
                                     .mnv_tag = mnv_tag,
                                     .report_mnv = ctx.report_mnv};
  AnnotatePhasing(merged, annotation, ordered_haps.empty() ? 0 : ordered_haps.front());

  // Key the merged record by the record index of its first component so it can be written in
  // position rather than appended at the end of the file.
  ctx.merged_by_first_idx.try_emplace(template_idx, std::move(merged));

  if (ctx.params.keep_components) {
    // Retain the original component records but mark them so downstream tools treat the merged
    // record as the canonical representation.
    for (const s32 vi : sorted_indices) {
      const s32 comp_idx = ctx.variants[static_cast<size_t>(vi)].vcf_record_index;
      ctx.out_records[static_cast<size_t>(comp_idx)]->SetFilter("PHASED_COMPONENT");
    }
  } else {
    for (const s32 vi : sorted_indices) {
      ctx.merged_away.insert(ctx.variants[static_cast<size_t>(vi)].vcf_record_index);
    }
  }

  Logging::Debug("Merged {} variants at {}:{} -> {}/{}",
                 sorted_indices.size(),
                 first_var.chrom,
                 first_var.pos + 1,
                 merged_ref,
                 merged_alt);
}

// Annotate a single phase block's records and merge it when configured to do so.
void ProcessPhaseBlock(WriteContext& ctx, const PhaseBlock& block) {
  if (block.variant_indices.size() < kMinBlockSize) {
    return;
  }

  vec<s32> sorted_indices;
  vec<s32> ordered_haps;
  OrderBlockByPosition(block, ctx.variants, sorted_indices, ordered_haps);

  const auto& first_var = ctx.variants[static_cast<size_t>(sorted_indices.front())];
  // Phase set (PS) id: the 1-based position of the leftmost variant in the block (VCF PS is 1-based).
  const s32 ps = static_cast<s32>(first_var.pos) + 1;
  const std::string pid = fmt::format("{}_{}", first_var.chrom, first_var.pos + 1);
  const std::string mnv_tag = ctx.report_mnv ? BuildMnvTag(sorted_indices, ctx.variants) : std::string{};
  const PhasingAnnotation annotation{.n_samples = ctx.n_samples,
                                     .tumor_idx = ctx.tumor_idx,
                                     .phase_set = ps,
                                     .pid = pid,
                                     .mnv_tag = mnv_tag,
                                     .report_mnv = ctx.report_mnv};

  // Annotate each constituent record with phasing fields.
  for (size_t k = 0; k < sorted_indices.size(); ++k) {
    const s32 vi = sorted_indices[k];
    const s32 rec_idx = ctx.variants[static_cast<size_t>(vi)].vcf_record_index;
    const s32 haplotype = ordered_haps.empty() ? 0 : ordered_haps[k];
    AnnotatePhasing(ctx.out_records[static_cast<size_t>(rec_idx)], annotation, haplotype);
  }

  const bool snv_only = IsSnvOnlyBlock(sorted_indices, ctx.variants);
  const bool do_merge = ctx.merge_enabled && (snv_only ? ctx.params.merge_mnv : ctx.params.merge_complex);
  if (do_merge) {
    MergeBlock(ctx, sorted_indices, ordered_haps, first_var, ps, pid, mnv_tag);
  }
}

// Clone input records against the phasing header before any annotations are added.
vec<io::VcfRecordPtr> CloneOutputRecords(const VcfData& vcf_data, const io::VcfHeaderPtr& hdr) {
  vec<io::VcfRecordPtr> out_records;
  out_records.reserve(vcf_data.records.size());
  for (const auto& record : vcf_data.records) {
    out_records.push_back(record->Clone(hdr));
  }
  return out_records;
}

// Open the BAM used to fill gaps between phased variants, if merging has been requested.
io::AlignmentReader OpenMergeBam(const Params& params, const bool merge_enabled) {
  io::AlignmentReader bam_handle;
  if (!merge_enabled || params.bam_path.empty()) {
    return bam_handle;
  }

  try {
    bam_handle = io::OpenAlignmentReader(params.bam_path);
  } catch (const std::exception& e) {  // NOSONAR(cpp:S1181)
    // A BAM is optional for merging; continue without read-derived gap filling when it cannot be opened.
    Logging::Warn("Failed to open BAM for merge gap-filling: {} ({})", params.bam_path.string(), e.what());
  }
  return bam_handle;
}

// Open the indexed FASTA used as a complete backbone for phased variant merging, if available.
std::unique_ptr<io::FastaReader> OpenMergeFasta(const Params& params, const bool merge_enabled) {
  if (!merge_enabled || params.ref_path.empty()) {
    return nullptr;
  }

  try {
    return std::make_unique<io::FastaReader>(params.ref_path.string());
  } catch (const std::exception& e) {  // NOSONAR(cpp:S1181) — reference failure degrades to BAM fallback
    Logging::Warn("Failed to open reference FASTA for merge gap-filling: {} ({})", params.ref_path.string(), e.what());
  }
  return nullptr;
}

// Annotate and optionally merge every phase block in the input.
void ProcessPhaseBlocks(WriteContext& ctx, const vec<PhaseBlock>& phase_blocks) {
  for (const auto& block : phase_blocks) {
    ProcessPhaseBlock(ctx, block);
  }
}

// Rescue filtered variants that were phased with a PASS variant.
void RescueFilteredVariants(const vec<PhaseBlock>& phase_blocks,
                            const vec<Variant>& variants,
                            vec<io::VcfRecordPtr>& out_records) {
  for (const auto& block : phase_blocks) {
    for (const s32 vi : block.variant_indices) {
      if (!variants[static_cast<size_t>(vi)].is_pass) {
        const s32 rec_idx = variants[static_cast<size_t>(vi)].vcf_record_index;
        out_records[static_cast<size_t>(rec_idx)]->SetFilter("PASS");
      }
    }
  }
}

// Write records in their original order and build the Tabix index after closing the BGZF stream.
void WriteRecords(const fs::path& output_path,
                  const io::VcfHeaderPtr& hdr,
                  const vec<io::VcfRecordPtr>& out_records,
                  const std::unordered_set<s32>& merged_away,
                  const std::unordered_map<s32, io::VcfRecordPtr>& merged_by_first_idx) {
  {
    const io::VcfWriter writer(output_path, hdr, io::NoIndex{});
    writer.WriteHeader();

    for (s32 rec_idx = 0; rec_idx < static_cast<s32>(out_records.size()); ++rec_idx) {
      // Emit a merged MNV record in the position of its first component.
      const auto merged_it = merged_by_first_idx.find(rec_idx);
      if (merged_it != merged_by_first_idx.end()) {
        writer.WriteRecord(merged_it->second);
      }
      if (merged_away.contains(rec_idx)) {
        continue;
      }
      writer.WriteRecord(out_records[static_cast<size_t>(rec_idx)]);
    }
  }

  if (output_path.extension() == ".gz" && tbx_index_build(output_path.c_str(), 0, &tbx_conf_vcf) != 0) {
    Logging::Warn("Failed to build Tabix index for '{}'", output_path.string());
  }
}

// Return the first non-missing value from a FORMAT string field. htslib concatenates values
// from multiple samples with NUL separators; the phase identifier is useful regardless of which
// sample column carries it.
std::string FirstPhaseValue(const std::string_view value) {
  size_t start = 0;
  while (start <= value.size()) {
    const size_t end = value.find('\0', start);
    std::string token{value.substr(start, end == std::string::npos ? std::string::npos : end - start)};
    if (!token.empty() && token != ".") {
      return token;
    }
    if (end == std::string::npos) {
      break;
    }
    start = end + 1;
  }
  return {};
}

// Read a pre-existing phase identifier from an input record. Mutect2 can provide PID without
// MNVTAG, while a previously phased VCF may only retain MNVTAG, so support both annotations.
std::string ExistingPhaseKey(const io::VcfRecordPtr& record) {
  if (record->HasFormatFieldNoCheck("PID")) {
    const std::string pid = FirstPhaseValue(record->GetFormatFieldStringNoCheck("PID"));
    if (!pid.empty()) {
      return "PID:" + pid;
    }
  }
  if (record->HasInfoFieldNoCheck("MNVTAG")) {
    const std::string mnv_tag = record->GetInfoFieldStringNoCheck("MNVTAG");
    if (!mnv_tag.empty() && mnv_tag != ".") {
      return "MNVTAG:" + mnv_tag;
    }
  }
  return {};
}

// Build phase blocks from phase identifiers already present in an input VCF. These annotations
// are emitted by Mutect2 and identify variants that were phased before this tool runs; retaining
// them lets the merger produce a complex record even when BAM evidence is unavailable.
vec<PhaseBlock> BuildExistingPhaseBlocks(const VcfData& vcf_data, const vec<Variant>& variants) {
  std::unordered_map<std::string, vec<s32>> indices_by_key;
  for (s32 vi = 0; vi < static_cast<s32>(variants.size()); ++vi) {
    const auto& variant = variants[static_cast<size_t>(vi)];
    const auto& record = vcf_data.records[static_cast<size_t>(variant.vcf_record_index)];
    const std::string key = ExistingPhaseKey(record);
    if (!key.empty()) {
      indices_by_key[key].push_back(vi);
    }
  }

  vec<PhaseBlock> blocks;
  for (auto& entry : indices_by_key) {
    auto& indices = entry.second;
    if (indices.size() < kMinBlockSize) {
      continue;
    }
    std::ranges::sort(indices, [&variants](const s32 lhs, const s32 rhs) {
      const auto& left = variants[static_cast<size_t>(lhs)];
      const auto& right = variants[static_cast<size_t>(rhs)];
      return left.pos < right.pos;
    });
    blocks.push_back(PhaseBlock{.variant_indices = std::move(indices)});
  }
  return blocks;
}

// Add input phase blocks that are not already represented by BAM-derived blocks. A pre-existing
// block may connect variants that were split across separate BAM-derived blocks, so it is skipped
// only when one supplied block already contains every one of its variants.
vec<PhaseBlock> AddExistingPhaseBlocks(const VcfData& vcf_data,
                                       const vec<PhaseBlock>& phase_blocks,
                                       const vec<Variant>& variants) {
  vec<PhaseBlock> result = phase_blocks;
  const auto existing_blocks = BuildExistingPhaseBlocks(vcf_data, variants);
  for (const auto& existing : existing_blocks) {
    const bool represented = std::ranges::any_of(result, [&existing](const PhaseBlock& block) {
      return std::ranges::all_of(existing.variant_indices, [&block](const s32 vi) {
        return std::ranges::find(block.variant_indices, vi) != block.variant_indices.end();
      });
    });
    if (!represented) {
      result.push_back(existing);
    }
  }
  return result;
}

}  // namespace

void WriteVcf(const fs::path& output_path,
              const VcfData& vcf_data,
              const vec<PhaseBlock>& phase_blocks,
              const vec<Variant>& variants,
              const Params& params) {
  const bool merge_enabled = params.merge_mnv || params.merge_complex;
  // Merging always emits an MNVTAG on the merged record, even when --report-mnv was not requested,
  // so that the merged representation records its constituent variants.
  const bool report_mnv = params.report_mnv || merge_enabled;

  const PhasingOutputOptions output_options{.merge_enabled = merge_enabled, .report_mnv = report_mnv};
  const io::VcfHeaderPtr hdr = BuildPhasedHeader(vcf_data, params, output_options);
  const vec<PhaseBlock> effective_phase_blocks =
      merge_enabled ? AddExistingPhaseBlocks(vcf_data, phase_blocks, variants) : phase_blocks;

  // Clone every input record against the new header so that appended FORMAT/INFO
  // definitions resolve when we set phasing fields.
  vec<io::VcfRecordPtr> out_records = CloneOutputRecords(vcf_data, hdr);

  // Records replaced by a merged MNV record (only when merging without keep-components).
  std::unordered_set<s32> merged_away;
  // Merged MNV records keyed by the original record index of their first component, so they can be
  // written in position (before that component's location) rather than appended at the end.
  std::unordered_map<s32, io::VcfRecordPtr> merged_by_first_idx;

  // Open optional resources used to reconstruct merged alleles.
  io::AlignmentReader bam_handle = OpenMergeBam(params, merge_enabled);

  // The reference provides a complete backbone across gaps that may not be covered by one read.
  std::unique_ptr<io::FastaReader> fasta_reader = OpenMergeFasta(params, merge_enabled);

  WriteContext ctx{.vcf_data = vcf_data,
                   .variants = variants,
                   .params = params,
                   .hdr = hdr,
                   .out_records = out_records,
                   .merged_away = merged_away,
                   .merged_by_first_idx = merged_by_first_idx,
                   .fasta_reader = fasta_reader.get(),
                   .bam_handle = bam_handle,
                   .n_samples = vcf_data.header->GetNumSamples(),
                   .tumor_idx = vcf_data.tumor_sample_idx,
                   .merge_enabled = merge_enabled,
                   .report_mnv = report_mnv};

  ProcessPhaseBlocks(ctx, effective_phase_blocks);

  // Rescue non-PASS variants that phased with a PASS variant.
  if (params.rescue_filtered) {
    RescueFilteredVariants(effective_phase_blocks, variants, out_records);
  }

  WriteRecords(output_path, hdr, out_records, merged_away, merged_by_first_idx);

  Logging::Info("Wrote phased VCF: {}", output_path.string());
}

}  // namespace xoos::svc::phaser
