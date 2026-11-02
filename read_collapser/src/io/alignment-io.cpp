#include "io/alignment-io.h"

#include <io/alignment.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>

#include "bam_cat.h"
#include "core/region.h"
#include "metrics/metrics.h"
#include "util/duplex-util.h"

namespace xoos::read_collapser {

bool HasRescuedSecondaryTag(const bam1_t* record) {
  const auto* tag_ptr = bam_aux_get(record, "YF");
  return tag_ptr != nullptr && bam_aux2i(tag_ptr) == 1;
}

// Determine if the alignment in @p record should be discarded based on @p options like mapping quality and alignment
// flags, increment appropriate values of @p metrics, and return true if the alignment should be discarded.
bool ShouldDiscardRecord(const ReadCollapserOptions& options, const Alignment* const alignment_ptr, Metrics& metrics) {
  bool discard = false;
  const auto* const record = alignment_ptr->record.get();
  ++metrics.input_records;
  if (record->core.qual < options.min_mapq) {
    ++metrics.discarded_low_mapq_records;
    discard = true;
  }
  if ((record->core.flag & options.exclude_flags) != 0) {
    // In markdup mode, secondary reads with YF:i:1 are rescued for clustering
    // instead of being discarded, but only if BAM_FSECONDARY is the sole
    // matching exclude flag. If any other excluded flag is set (e.g., QC fail),
    // the read is still discarded.
    const auto matched_flags = record->core.flag & options.exclude_flags;
    const bool matched_only_secondary = matched_flags == BAM_FSECONDARY;
    const bool is_rescued =
        options.markdup_options.cluster_rescued_secondaries && matched_only_secondary && HasRescuedSecondaryTag(record);
    if (!is_rescued) {
      ++metrics.discarded_by_flags_records;
      discard = true;
    }
  }
  if (options.max_discordant_duplex_error_percentage) {
    const auto duplex_error_rate = GetDiscordantDuplexErrorRate(record);
    if (duplex_error_rate && *duplex_error_rate * 100.0 > *options.max_discordant_duplex_error_percentage) {
      ++metrics.discarded_high_discordant_duplex_percentage_records;
      discard = true;
    }
  }
  if (options.cluster_by_umi) {
    const bool has_umi5p = alignment_ptr->umi5p.has_value();
    const bool has_umi3p = alignment_ptr->umi3p.has_value();
    if (!has_umi5p && !has_umi3p) {
      ++metrics.discarded_missing_umi_records;
      discard = true;
    }
  }
  if (discard) {
    ++metrics.discarded_total_records;
  }
  return discard;
}

void ReadAlignments(const ReadCollapserOptions& options,
                    const io::AlignmentReader& reader,
                    const Region& region,
                    vec<AlignmentPtr>& alignments,
                    Metrics& metrics) {
  const auto itr = io::SamItrQueryI(reader.idx.get(), region.tid, region.start, region.end);
  while (true) {
    auto record = io::Bam1Ptr{bam_init1()};
    const auto has_more_reads = io::SamItrNext(reader.bam_file.get(), itr.get(), record.get());
    if (!has_more_reads) {
      break;
    }
    // Skip alignments with a start position outside of the region. This is to avoid double counting
    // reads which start in one region but overlap another region.
    // TODO: This should be improved to handle alignments which start outside of the region,
    //  but overlap the region. As it is now these alignments will not be included. This will have a
    //  very small affect on some target panels, where some alignments that start outside of the region
    //  but within the wiggle room are not included.
    if (record->core.pos < region.start || record->core.pos >= region.end) {
      continue;
    }
    const auto alignment_ptr =
        std::make_shared<Alignment>(std::move(record), options.cluster_by_umi, options.ignore_read_name_parsing_errors);
    if (ShouldDiscardRecord(options, alignment_ptr.get(), metrics)) {
      continue;
    }
    alignments.emplace_back(alignment_ptr);
  }
}

HtsIdxWriter::HtsIdxWriter(HtsIdxWriter&& other) noexcept : bam(other.bam), idx_fn(std::move(other.idx_fn)) {
  // Nullify the other bam pointer to prevent double deletion.
  other.bam = nullptr;
}

HtsIdxWriter& HtsIdxWriter::operator=(HtsIdxWriter&& other) noexcept {
  if (this != &other) {
    // First, clean up any existing resources.
    if (bam != nullptr) {
      io::HtsFlush(bam);
      io::SamIdxSave(bam);
    }
    // Move the resources from other to this.
    bam = other.bam;
    idx_fn = std::move(other.idx_fn);
    // Nullify the other bam pointer to prevent double deletion.
    other.bam = nullptr;
  }
  return *this;
}

HtsIdxWriter::~HtsIdxWriter() {
  // Add a flush first to update the timestamp of the BAM if necessary,
  // this helps to avoid a samtools warning about the index being older than the BAM.
  if (bam != nullptr) {
    io::HtsFlush(bam);
    io::SamIdxSave(bam);
  }
}

AlignmentWriter OpenAlignmentWriter(const fs::path& location, sam_hdr_t* hdr, const bool write_index) {
  htsFormat fmt;
  fmt.category = sequence_data;
  fmt.format = bam;
  fmt.compression = bgzf;
  fmt.compression_level = 1;
  fmt.specific = nullptr;

  create_directories(location.parent_path());
  auto bam = io::HtsOpenFormat(location, "wbz1", &fmt);
  std::optional<HtsIdxWriter> bam_idx;

  if (write_index) {
    bam_idx.emplace(bam.get(), location.string() + ".bai");
    io::SamIdxInit(bam.get(), hdr, 0, bam_idx.value().idx_fn);
  }

  io::SamHdrWrite(bam.get(), hdr);

  return AlignmentWriter{.bam = std::move(bam), .bam_idx = std::move(bam_idx)};
}

void RemovePgTag(bam1_t* const record) {
  auto* const pg_tag = bam_aux_get(record, "PG");
  if (pg_tag != nullptr) {
    bam_aux_del(record, pg_tag);
  }
}

void WriteAlignment(const bam1_t* const record, const bool mark_duplicate, const io::HtsFilePtr& output) {
  WriteAlignment(record, mark_duplicate, output, false);
}

void WriteAlignment(const bam1_t* const record,
                    const bool mark_duplicate,
                    const io::HtsFilePtr& output,
                    const bool remove_pg_tag) {
  auto record_modified = io::Bam1Ptr{bam_dup1(record)};
  // Clear the existing duplicate flag if it exists, since we will set it based on our clustering results.
  // This is to avoid the case where the input BAM has been duplicate marked by another tool.
  record_modified->core.flag &= static_cast<u16>(~static_cast<u32>(BAM_FDUP));
  if (mark_duplicate) {
    record_modified->core.flag |= BAM_FDUP;
  }
  if (remove_pg_tag) {
    RemovePgTag(record_modified.get());
  }
  io::SamWrite1(output.get(), output->bam_header, record_modified.get());
}

void WriteAlignment(const bam1_t* const record,
                    const bool mark_duplicate,
                    const ClusterId& cluster_id,
                    const u32 cluster_size,
                    const io::HtsFilePtr& output) {
  WriteAlignment(record, mark_duplicate, cluster_id, cluster_size, output, false);
}

void WriteAlignment(const bam1_t* const record,
                    const bool mark_duplicate,
                    const ClusterId& cluster_id,
                    const u32 cluster_size,
                    const io::HtsFilePtr& output,
                    const bool remove_pg_tag) {
  auto record_modified = io::Bam1Ptr{bam_dup1(record)};
  // Clear the existing duplicate flag if it exists, since we will set it based on our clustering results.
  // This is to avoid the case where the input BAM has been duplicate marked by another tool.
  record_modified->core.flag &= static_cast<u16>(~static_cast<u32>(BAM_FDUP));
  if (mark_duplicate) {
    record_modified->core.flag |= BAM_FDUP;
  }
  if (remove_pg_tag) {
    RemovePgTag(record_modified.get());
  }
  io::BamAuxAppend(record_modified.get(), "DI", fmt::format("{}:{}", cluster_id.super_region_id, cluster_id.count));
  io::BamAuxAppend(record_modified.get(), "DS", cluster_size);
  io::SamWrite1(output.get(), output->bam_header, record_modified.get());
}

void WriteAlignments(const vec<AlignmentPtr>& alignments,
                     const bool remove_duplicates,
                     const bool include_cluster_info,
                     const io::HtsFilePtr& output) {
  WriteAlignments(alignments, remove_duplicates, include_cluster_info, output, false);
}

void WriteAlignments(const vec<AlignmentPtr>& alignments,
                     const bool remove_duplicates,
                     const bool include_cluster_info,
                     const io::HtsFilePtr& output,
                     const bool remove_pg_tag) {
  for (const auto& alignment : alignments) {
    if (remove_duplicates && alignment->is_duplicate) {
      continue;
    }
    const auto mark_duplicate = !remove_duplicates && alignment->is_duplicate;

    if (auto* cluster = alignment->cluster; include_cluster_info && cluster != nullptr) {
      auto cluster_size = static_cast<u32>(cluster->alignments.size());
      WriteAlignment(alignment->record.get(), mark_duplicate, cluster->cluster_id, cluster_size, output, remove_pg_tag);
    } else {
      WriteAlignment(alignment->record.get(), mark_duplicate, output, remove_pg_tag);
    }
  }
}

void ReadAndWriteUnmappedReads(const ReadCollapserOptions& options,
                               const CreateClusterId& create_cluster_id,
                               const io::AlignmentReader& reader,
                               const UnmappedAlignmentWriter& writer,
                               Metrics& metrics) {
  const auto itr = io::HtsItrPtr{sam_itr_queryi(reader.idx.get(), HTS_IDX_NOCOOR, 0, 0)};
  if (itr == nullptr) {
    // If the region is unmapped and the bam file does not contain any reads, itr will be nullptr.
    // Otherwise, as long as the bam file contains at least one read (mapped or unmapped), itr will be a valid pointer
    // for an unmapped region. This is a limitation of sam_itr_queryi, which has special handling of HTS_IDX_NOCOOR and
    // assumes all unmapped reads can be found at the end of the bam file (i.e. when a bam file is empty, this section
    // cannot exist).
    return;
  }
  while (true) {
    auto record = io::Bam1Ptr{bam_init1()};
    const auto has_more_reads = io::SamItrNext(reader.bam_file.get(), itr.get(), record.get());
    if (!has_more_reads) {
      break;
    }
    const auto alignment_ptr =
        std::make_shared<Alignment>(std::move(record), options.cluster_by_umi, options.ignore_read_name_parsing_errors);
    if (ShouldDiscardRecord(options, alignment_ptr.get(), metrics)) {
      continue;
    }
    ++metrics.unmapped_reads;
    writer(create_cluster_id(), alignment_ptr->record.get());
  }
}

void ConcatenateBamFiles(const vec<fs::path>& input_bam_files, const fs::path& output_bam_file) {
  const char* first_file = input_bam_files[0].c_str();
  vec<const char*> input_filenames;
  for (const auto& input_bam_file : input_bam_files) {
    input_filenames.emplace_back(input_bam_file.c_str());
  }
  const auto file_count = static_cast<s32>(input_bam_files.size());
  // bam_cat will take ownership of first_sam_file and header
  samFile* first_sam_file = io::HtsOpen(first_file, "r").release();
  const char* output_file = output_bam_file.c_str();
  const s32 bam_cat_result =
      bam_cat(first_sam_file, file_count, input_filenames.data(), nullptr, output_file, nullptr, 1);
  if (bam_cat_result < 0) {
    throw error::Error("bam_cat failed");
  }
}

}  // namespace xoos::read_collapser
