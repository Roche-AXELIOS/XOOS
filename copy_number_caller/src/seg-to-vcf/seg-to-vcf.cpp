#include "seg-to-vcf/seg-to-vcf.h"

#include <algorithm>
#include <cmath>
#include <fstream>

#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <xoos/enum/enum-util.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/metadata-util.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>
#include <xoos/util/missing-value.h>

#include "cli/copy-number-caller-options.h"
#include "fmt/core.h"
#include "io/copy-number-caller-default-filenames.h"
#include "likelihood/likelihood-flags.h"
#include "segmentation/read-segments.h"
#include "segmentation/segment-type.h"
#include "segmentation/segments-header.h"

namespace xoos::cnc {

// The "Main" function only supports germline for now
void WriteSegmentsToVcfMain(const cnc::CopyNumberCallerOptions& options) {
  const auto vcf_out_fname = options.output_dir / kDefaultGermlineCNCallsetVcfOutput;
  file::CheckOutputFilePathsForWriting({vcf_out_fname}, options.overwrite);
  WriteSegmentsToVcfParams params;
  if (!options.reference_genome_fai_fname.empty()) {
    params.seq_lengths = ParseFai(options.reference_genome_fai_fname);
    params.seq_order = GetContigOrder(options.reference_genome_fai_fname);
  }
  params.sample_id = options.sample_metadata_options.sample_id;
  params.sex = options.sample_metadata_options.sex.value();
  const vec<GenomicSegment> segments(ReadSegments(options.segments_fname.value(), SegmentType::kGermlineLikelihood));
  const SegmentsHeader segments_header = ReadHeaderFromSegments(options.segments_fname.value());
  params.reference_file =
      segments_header.reference_file.empty() ? options.reference_genome_fname.string() : segments_header.reference_file;
  params.mode = CopyNumberCallerModes::kGermlineNormalWGS;
  params.command_line_info = options.command_line_info;
  WriteSegmentsToVcf(vcf_out_fname, segments, params);
}

/**
 * @brief Determines if a genomic segment is located on the male allosome (sex chromosome)
 *        and is not within the pseudoautosomal region (PAR).
 *
 * @param seg The genomic segment to evaluate.
 * @return true if the segment is on the allosome, belongs to a male, and is not in the PAR; false otherwise.
 */
namespace {
bool IsMaleAllosomeNotInPar(const GenomicSegment& seg) {
  return seg.in_allosome.value() && seg.sex == sex_predict::Sex::kMale && !seg.in_pseudo_autosomal_region.value();
}
}  // namespace

/**
 * @brief Determines the copy number status for a male allosome segment not in the pseudoautosomal region (PAR).
 *
 * This function evaluates the total copy number of a given genomic segment and returns the corresponding
 * CopyNumberStatus:
 *   - Returns CopyNumberStatus::kRef if the total copy number is 1 (reference state for male allosomes).
 *   - Returns CopyNumberStatus::kGain if the total copy number is greater than 1 (copy number gain).
 *   - Returns CopyNumberStatus::kLoss if the total copy number is less than 1 (copy number loss).
 *
 * @param seg The GenomicSegment representing the region to evaluate.
 * @return CopyNumberStatus The determined copy number status for the segment.
 */
namespace {
CopyNumberStatus GetMaleAllosomeNotInParCopyNumberStatus(const GenomicSegment& seg) {
  if (seg.total_copy_number == 1) {
    return CopyNumberStatus::kRef;
  }
  if (seg.total_copy_number > 1) {
    return CopyNumberStatus::kGain;
  }
  return CopyNumberStatus::kLoss;
}
}  // namespace

/**
 * @brief Determines the copy number status of a given genomic segment.
 *
 * This function evaluates the copy number status for the provided GenomicSegment.
 * If the segment is a male allosome not in the pseudoautosomal region (PAR),
 * it delegates the status determination to GetMaleAllosomeNotInParCopyNumberStatus.
 * Otherwise, it classifies the segment as reference (kRef) if the total copy number is 2,
 * as gain (kGain) if the total copy number is greater than 2, or as loss (kLoss) if less than 2.
 *
 * @param seg The genomic segment to evaluate.
 * @return CopyNumberStatus The determined copy number status for the segment.
 */
CopyNumberStatus GetCopyNumberStatus(const GenomicSegment& seg) {
  using enum CopyNumberStatus;
  if (IsMaleAllosomeNotInPar(seg)) {
    return GetMaleAllosomeNotInParCopyNumberStatus(seg);
  }
  if (seg.total_copy_number == kDiploidCopyNumber) {
    return CopyNumberStatus::kRef;
  }
  if (seg.total_copy_number > kDiploidCopyNumber) {
    return CopyNumberStatus::kGain;
  }
  return CopyNumberStatus::kLoss;
}

/**
 * @brief Determines the somatic copy number status for a given genomic segment.
 *
 * This function evaluates the copy number status of the provided GenomicSegment.
 * For male allosomes not in the pseudoautosomal region (PAR), it delegates to a specialized handler.
 * For other cases, it classifies the segment as loss, reference, copy-neutral loss of heterozygosity (cnLOH),
 * gain, or gain with LOH, based on the total and minor copy numbers.
 *
 * @param seg The GenomicSegment to evaluate.
 * @return CopyNumberStatus The determined somatic copy number status.
 * @throws std::runtime_error If the segment does not match any valid copy number status.
 */
CopyNumberStatus GetSomaticCopyNumberStatus(const GenomicSegment& seg) {
  using enum CopyNumberStatus;
  if (IsMaleAllosomeNotInPar(seg)) {
    return GetMaleAllosomeNotInParCopyNumberStatus(seg);
  }
  if (seg.total_copy_number <= 1) {
    return CopyNumberStatus::kLoss;
  }
  if (seg.total_copy_number == kDiploidCopyNumber) {
    if (seg.minor_copy_number.has_value()) {
      if (seg.minor_copy_number.value() == 0) {
        return CopyNumberStatus::kCnLoh;
      }
    }
    return CopyNumberStatus::kRef;
  }
  if (seg.total_copy_number > kDiploidCopyNumber) {
    if (seg.minor_copy_number.has_value()) {
      if (seg.minor_copy_number.value() == 0) {
        return CopyNumberStatus::kGainLoh;
      }
    }
    return CopyNumberStatus::kGain;
  }
  throw std::runtime_error("invalid somatic copy number status");
}

namespace {
std::string ConstructVcfIdFieldFromGenomicSegment(const GenomicSegment& seg, CopyNumberStatus copy_number_status) {
  using enum CopyNumberStatus;
  std::string ret = "XOOS";
  switch (copy_number_status) {
    case kRef:
      ret += ":REF";
      break;
    case kGain:
      ret += ":GAIN";
      break;
    case kLoss:
      ret += ":LOSS";
      break;
    case kCnLoh:
      ret += ":CNLOH";
      break;
    case kGainLoh:
      ret += ":GAINLOH";
      break;
    default:
      throw std::runtime_error("invalid copy number status");
  }
  ret += fmt::format(":{}:{}-{}", seg.contig, seg.start + 1, seg.end);
  return ret;
}

std::string ConstructVcfRefFieldFromGenomicSegment() {
  return "N";
}

std::string ConstructVcfAltFieldFromGenomicSegment(CopyNumberStatus copy_number_status) {
  using enum CopyNumberStatus;
  switch (copy_number_status) {
    case kRef:
      return std::string{util::kMissingSpecField};
    case kGain:
      return "<DUP>";
    case kLoss:
      return "<DEL>";
    case kCnLoh:
    case kGainLoh:
      return "<DEL>,<DUP>";
    default:
      throw std::runtime_error("invalid copy number status");
  }
}

std::string ConstructVcfQualFieldFromGenomicSegment() {
  return std::string{util::kMissingSpecField};
}
}  // namespace

std::string ConstructVcfFilterFieldFromGenomicSegment(const GenomicSegment& seg) {
  if (!seg.flags.has_value() || seg.flags.value().empty()) {
    return kLikelihoodFlagPass;
  }
  return FlagsToString(seg.flags.value());
}

namespace {
std::string ConstructVcfInfoFieldFromGenomicSegment(const GenomicSegment& seg,
                                                    const CopyNumberStatus copy_number_status,
                                                    const CopyNumberCallerModes mode) {
  using enum CopyNumberStatus;
  std::string ret;
  size_t end = seg.end;
  size_t reflen = seg.end - seg.start;
  switch (copy_number_status) {
    case kRef:
      ret += "";
      break;
    case kGain:
      ret += fmt::format("SVLEN={};SVTYPE=CNV;", reflen);
      break;
    case kGainLoh:
      ret += fmt::format("SVLEN=-{},{};SVTYPE=CNV;", reflen, reflen);
      break;
    case kLoss:
      ret += fmt::format("SVLEN=-{};SVTYPE=CNV;", reflen);
      break;
    case kCnLoh:
      ret += fmt::format("SVLEN=-{},{};SVTYPE=CNV;", reflen, reflen);
      break;
    default:
      throw std::runtime_error("invalid copy number status");
  }
  ret += fmt::format("END={};REFLEN={}", end, reflen);
  // Only tumor-only mode classifies segments as germline CNVs and declares the GERMLINE INFO header,
  // so the flag must only be emitted here to avoid an undeclared INFO field in other modes.
  if ((mode == CopyNumberCallerModes::kSomaticTumorOnlyWGS) && seg.is_germline_cnv.value_or(false)) {
    ret += ";GERMLINE";
  }
  return ret;
}
}  // namespace

/**
 * @brief Constructs a VCF FORMAT field string from a genomic segment with copy number information.
 *
 * This function generates a VCF FORMAT field that includes genotype and copy number statistics
 * based on the provided genomic segment data and copy number status. The output format varies
 * depending on the copy number caller mode (germline vs somatic analysis).
 *
 * @param seg The genomic segment containing copy number and statistical data
 * @param copy_number_status The copy number status (reference, gain, loss, CN-LOH, etc.)
 * @param mode The copy number caller mode (germline normal WGS or somatic tumor-normal WGS)
 *
 * @return A formatted string containing the VCF FORMAT field header and corresponding values
 *         separated by tabs. For germline mode, includes GT:MRDR:LR:NUMTARGETS:MMMAPQ:CN.
 *         For somatic mode, includes NUMTARGETS:MMMAPQ:MRDR:LR:NUMSNPS:MBAF:CN:MCN:GT.
 *
 * @throws std::runtime_error If copy_number_status is invalid or mode is unsupported
 *
 * @note Genotype assignment logic:
 *       - Reference: "0/0"
 *       - Gain: "./1", "0/1", or "1/1" depending on minor copy number availability and value
 *       - Loss: "1" (allosome), "1/1" (homozygous deletion), or "0/1" (heterozygous deletion)
 *       - CN-LOH/Gain-LOH: "1/2" (representing <DEL>/<DUP>)
 */
namespace {
// Determine the VCF genotype string for a segment given its copy-number status.
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
std::string ComputeGenotypeString(const GenomicSegment& seg, const CopyNumberStatus copy_number_status) {
  switch (copy_number_status) {
    case CopyNumberStatus::kRef:
      return "0/0";
    case CopyNumberStatus::kGain:
      // If we lack minor copy number data, then we cannot indicate which allele was gained
      if (!seg.minor_copy_number.has_value()) {
        return "./1";
      }
      // If minor copy number is 1, then we know only one allele was duplicated. If it is 2, then we know
      // both alleles were duplicated
      if (seg.minor_copy_number.value() == 1) {
        return "0/1";
      }
      if (seg.minor_copy_number.value() >= 2) {
        return "1/1";
      }
    case CopyNumberStatus::kLoss:
      if (seg.in_allosome.value() && !seg.in_pseudo_autosomal_region.value()) {
        return "1";
      }
      if (seg.total_copy_number == 0) {
        // Homozygous deletion
        return "1/1";
      }
      // Heterozygous deletion
      return "0/1";
    case CopyNumberStatus::kCnLoh:
    case CopyNumberStatus::kGainLoh:
      // First allele will be <DEL> and second allele will be <DUP>
      return "1/2";
    default:
      throw std::runtime_error("invalid copy number status");
  }
}

// Add the somatic-specific FORMAT fields (MCN, NUMSNPS, MBAF) to the shared field map.
void AddSomaticFormatFields(const GenomicSegment& seg, StrMap<std::string>& format_fields) {
  std::string minor_cn{util::kMissingSpecField};
  if (seg.minor_copy_number.has_value()) {
    minor_cn = std::to_string(seg.minor_copy_number.value());
  }

  size_t num_snps = 0;
  if (seg.num_snps.has_value()) {
    num_snps = seg.num_snps.value();
  }

  std::string mbaf{util::kMissingSpecField};
  if (seg.mbaf.has_value()) {
    mbaf = fmt::format("{:.2f}", seg.mbaf.value());
  }

  format_fields["MCN"] = minor_cn;
  format_fields["NUMSNPS"] = std::to_string(num_snps);
  format_fields["MBAF"] = mbaf;
}

std::string ConstructVcfFormatFieldFromGenomicSegment(const GenomicSegment& seg,
                                                      const CopyNumberStatus copy_number_status,
                                                      const CopyNumberCallerModes mode) {
  const std::string gt = ComputeGenotypeString(seg, copy_number_status);
  f64 mrdr = std::pow(2, seg.mean_logr.value());
  f64 lr = seg.mean_logr.value();
  const size_t cn = seg.total_copy_number.value();
  std::string mmmapq{util::kMissingSpecField};
  if (seg.avg_mean_mapq.has_value()) {
    mmmapq = fmt::format("{:.2f}", seg.avg_mean_mapq.value());
  }

  // Create map of format fields shared between different modes.
  // We will need more to this depending on the mode
  StrMap<std::string> format_fields = {{"GT", gt},
                                       {"MRDR", fmt::format("{:.2f}", mrdr)},
                                       {"LR", fmt::format("{:.2f}", lr)},
                                       {"NUMTARGETS", std::to_string(seg.num_obs.value())},
                                       {"MMMAPQ", mmmapq},
                                       {"CN", std::to_string(cn)}};

  // Add mode-specific fields
  if (IsSomaticLikeMode(mode)) {
    AddSomaticFormatFields(seg, format_fields);
  }

  // Define format field order based on mode
  // GT must be the first subfield
  vec<std::string> field_order;
  if (mode == CopyNumberCallerModes::kGermlineNormalWGS) {
    field_order = {"GT", "NUMTARGETS", "MMMAPQ", "MRDR", "LR", "CN"};
  } else if (IsSomaticLikeMode(mode)) {
    field_order = {"GT", "NUMTARGETS", "MMMAPQ", "MRDR", "LR", "NUMSNPS", "MBAF", "CN", "MCN"};
  } else {
    throw std::runtime_error("unsupported copy number caller mode for VCF FORMAT field construction");
  }

  // Build format line
  std::string ret;
  std::string values = "\t";
  for (size_t i = 0; i < field_order.size(); ++i) {
    const std::string& field = field_order[i];
    ret += field;
    values += format_fields[field];
    if (i < field_order.size() - 1) {
      ret += ":";
      values += ":";
    }
  }
  ret += values;

  return ret;
}

std::string SegmentToRecordString(const GenomicSegment& seg, const CopyNumberCallerModes mode) {
  CopyNumberStatus copy_number_status = CopyNumberStatus::kUnknown;
  if (IsSomaticLikeMode(mode)) {
    copy_number_status = GetSomaticCopyNumberStatus(seg);
  } else {
    copy_number_status = GetCopyNumberStatus(seg);
  }
  std::string id = ConstructVcfIdFieldFromGenomicSegment(seg, copy_number_status);
  std::string ref = ConstructVcfRefFieldFromGenomicSegment();
  std::string alt = ConstructVcfAltFieldFromGenomicSegment(copy_number_status);
  std::string qual = ConstructVcfQualFieldFromGenomicSegment();
  std::string filter{util::kMissingSpecField};
  // not currently outputting filters for Somatic at the moment, but we can if requested, since the flags field should
  // already be populated
  if (mode == CopyNumberCallerModes::kGermlineNormalWGS) {
    filter = ConstructVcfFilterFieldFromGenomicSegment(seg);
  }
  std::string info = ConstructVcfInfoFieldFromGenomicSegment(seg, copy_number_status, mode);
  std::string format = ConstructVcfFormatFieldFromGenomicSegment(seg, copy_number_status, mode);
  return fmt::format(
      "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", seg.contig, seg.start + 1, id, ref, alt, qual, filter, info, format);
}

// Helper functions to get metadata lines by ID
std::string GetAltMetadataLine(const std::string& id) {
  if (id == "CNV") {
    return "##ALT=<ID=CNV,Description=\"Copy number variant region\">";
  }
  if (id == "DEL") {
    return "##ALT=<ID=DEL,Description=\"Region of lowered copy number relative to the reference, or a deletion "
           "breakpoint\">";
  }
  if (id == "DUP") {
    return "##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the reference, or a tandem "
           "duplication breakpoint\">";
  }
  throw std::runtime_error(fmt::format("Unknown ALT ID: {}", id));
}

std::string GetInfoMetadataLine(const std::string& id) {
  if (id == "REFLEN") {
    return "##INFO=<ID=REFLEN,Number=1,Type=Integer,Description=\"Number of REF positions included in this record\">";
  }
  if (id == "SVLEN") {
    return "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">";
  }
  if (id == "SVTYPE") {
    return "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">";
  }
  if (id == "END") {
    return "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">";
  }
  if (id == "GERMLINE") {
    return "##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description=\"Segment fits the germline copy-number model better "
           "than the somatic model (tumor-only mode)\">";
  }
  throw std::runtime_error(fmt::format("Unknown INFO ID: {}", id));
}

std::string GetFilterMetadataLine(const std::string& id) {
  if (id == "lowMmMapq") {
    return "##FILTER=<ID=lowMmMapq,Description=\"Low mean of mean MAPQ for overlapping targets\">";
  }
  if (id == "cnvLength") {
    return "##FILTER=<ID=cnvLength,Description=\"CNV with length below threshold\">";
  }
  if (id == "nonIntegerTcn") {
    return "##FILTER=<ID=nonIntegerTcn,Description=\"Copy-number call with observed LR that deviates too far from "
           "the expected LR of an integer total copy number value\">";
  }
  throw std::runtime_error(fmt::format("Unknown FILTER ID: {}", id));
}
}  // namespace

/**
 * @brief Returns the VCF FORMAT metadata line corresponding to the given FORMAT ID.
 *
 * This function generates a VCF FORMAT metadata line string for a specified FORMAT field ID.
 * Supported IDs include "GT", "MRDR", "LR", "NUMTARGETS", "NUMSNPS", "MMMAPQ", "CN", and "MCN", each corresponding
 * to a specific FORMAT field in the VCF specification. If an unknown ID is provided, the function
 * throws a std::runtime_error.
 *
 * @param id The FORMAT field ID for which to generate the metadata line.
 * @return std::string The corresponding VCF FORMAT metadata line.
 * @throws std::runtime_error If the provided ID is not recognized.
 */
namespace {
std::string GetFormatMetadataLine(const std::string& id) {
  if (id == "GT") {
    return "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
  }
  if (id == "MRDR") {
    return "##FORMAT=<ID=MRDR,Number=1,Type=Float,Description=\"Mean of the read depth ratios of overlapping "
           "targets\">";
  }
  if (id == "LR") {
    return "##FORMAT=<ID=LR,Number=1,Type=Float,Description=\"Log2 of the MRDR\">";
  }
  if (id == "NUMTARGETS") {
    return "##FORMAT=<ID=NUMTARGETS,Number=1,Type=Integer,Description=\"Number of overlapping targets\">";
  }
  if (id == "MMMAPQ") {
    return "##FORMAT=<ID=MMMAPQ,Number=1,Type=Float,Description=\"Mean of mean MAPQ of overlapping targets\">";
  }
  if (id == "CN") {
    return "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Estimated total copy number\">";
  }
  if (id == "MCN") {
    return "##FORMAT=<ID=MCN,Number=1,Type=Integer,Description=\"Estimated minor copy number\">";
  }
  if (id == "NUMSNPS") {
    return "##FORMAT=<ID=NUMSNPS,Number=1,Type=Integer,Description=\"Number of overlapping SNPs\">";
  }
  if (id == "MBAF") {
    return "##FORMAT=<ID=MBAF,Number=1,Type=Float,Description=\"Mean mirrored B-allele fraction of overlapping SNPs. "
           "SNPs with BAF > 0.5 are mirrored.\">";
  }
  throw std::runtime_error(fmt::format("Unknown FORMAT ID: {}", id));
}

std::string GetVcfHeaderLine(const std::string& sample_id) {
  return fmt::format("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", sample_id);
}

// Write the sample-level header lines (purity, ploidy, and related warnings) for somatic-like modes.
template <typename LineWriter>
void WriteSomaticSampleHeaders(const WriteSegmentsToVcfParams& params, const LineWriter& write_line) {
  if (!params.purity.has_value() || !params.ploidy.has_value()) {
    throw std::runtime_error("Purity and ploidy must be provided for somatic tumor-normal and tumor-only WGS modes");
  }
  write_line(fmt::format("##Purity={:.2f}", params.purity.value()));
  write_line(fmt::format("##Ploidy={:.2f}", params.ploidy.value()));
  if (params.vcf_purity_source.has_value()) {
    write_line(fmt::format("##PuritySource={}", enum_util::FormatEnumName(params.vcf_purity_source.value())));
  }
  if (params.total_copy_number_prior.has_value()) {
    write_line(
        fmt::format("##TotalCopyNumberPrior={}", enum_util::FormatEnumName(params.total_copy_number_prior.value())));
  }
  if (params.has_high_extreme_baf_proportion) {
    write_line(std::string(kSampleWarningHighExtremeBafProportion));
  }
}

// Write contig header lines ordered by their sequence index.
template <typename LineWriter>
void WriteContigHeaders(const WriteSegmentsToVcfParams& params, const LineWriter& write_line) {
  vec<std::tuple<size_t, std::string>> contig_header_vec;
  for (const auto& [contig, index] : params.seq_order) {
    if (params.seq_lengths.contains(contig)) {
      size_t length = params.seq_lengths.at(contig);
      contig_header_vec.emplace_back(index, fmt::format("##contig=<ID={},length={}>", contig, length));
    }
  }
  std::ranges::sort(contig_header_vec);
  for (const auto& [index, contig_header_line] : contig_header_vec) {
    write_line(contig_header_line);
  }
}

// Write the ALT/FILTER/INFO/FORMAT metadata lines for the germline normal WGS mode.
template <typename LineWriter>
void WriteGermlineMetadataLines(const LineWriter& write_line) {
  write_line(GetAltMetadataLine("CNV"));
  write_line(GetAltMetadataLine("DEL"));
  write_line(GetAltMetadataLine("DUP"));
  write_line(GetFilterMetadataLine("lowMmMapq"));
  write_line(GetFilterMetadataLine("cnvLength"));
  write_line(GetFilterMetadataLine("nonIntegerTcn"));
  write_line(GetInfoMetadataLine("REFLEN"));
  write_line(GetInfoMetadataLine("SVLEN"));
  write_line(GetInfoMetadataLine("SVTYPE"));
  write_line(GetInfoMetadataLine("END"));
  write_line(GetFormatMetadataLine("GT"));
  write_line(GetFormatMetadataLine("MRDR"));
  write_line(GetFormatMetadataLine("LR"));
  write_line(GetFormatMetadataLine("NUMTARGETS"));
  write_line(GetFormatMetadataLine("MMMAPQ"));
  write_line(GetFormatMetadataLine("CN"));
}

// Write the ALT/FILTER/INFO/FORMAT metadata lines for somatic-like modes.
template <typename LineWriter>
void WriteSomaticMetadataLines(const CopyNumberCallerModes mode, const LineWriter& write_line) {
  write_line(GetAltMetadataLine("CNV"));
  write_line(GetAltMetadataLine("DEL"));
  write_line(GetAltMetadataLine("DUP"));
  write_line(GetFilterMetadataLine("cnvLength"));
  write_line(GetInfoMetadataLine("REFLEN"));
  write_line(GetInfoMetadataLine("SVLEN"));
  write_line(GetInfoMetadataLine("SVTYPE"));
  write_line(GetInfoMetadataLine("END"));
  // Tumor-only mode flags segments that fit the germline copy-number model better than the somatic model.
  if (mode == CopyNumberCallerModes::kSomaticTumorOnlyWGS) {
    write_line(GetInfoMetadataLine("GERMLINE"));
  }
  write_line(GetFormatMetadataLine("GT"));
  write_line(GetFormatMetadataLine("NUMTARGETS"));
  write_line(GetFormatMetadataLine("MMMAPQ"));
  write_line(GetFormatMetadataLine("MRDR"));
  write_line(GetFormatMetadataLine("LR"));
  write_line(GetFormatMetadataLine("NUMSNPS"));
  write_line(GetFormatMetadataLine("MBAF"));
  write_line(GetFormatMetadataLine("CN"));
  write_line(GetFormatMetadataLine("MCN"));
}

// Write all VCF metadata (##) header lines for the given mode.
template <typename LineWriter>
void WriteVcfMetadataLines(const WriteSegmentsToVcfParams& params, const LineWriter& write_line) {
  if (params.mode == CopyNumberCallerModes::kGermlineNormalWGS) {
    WriteGermlineMetadataLines(write_line);
  } else if (IsSomaticLikeMode(params.mode)) {
    WriteSomaticMetadataLines(params.mode, write_line);
  } else {
    throw std::runtime_error("Unsupported CopyNumberCaller mode for writing segments to VCF");
  }
}
}  // namespace

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void WriteSegmentsToVcf(const fs::path& fname,
                        const vec<GenomicSegment>& segments,
                        const WriteSegmentsToVcfParams& params) {
  // Check if file should be bgzipped based on extension
  const std::string extension = fname.extension().string();
  const bool use_bgzip = (extension == ".gz" || extension == ".bgz");

  io::BgzfPtr bgzf_fp;
  std::ofstream ofs;

  if (use_bgzip) {
    bgzf_fp.reset(bgzf_open(fname.c_str(), "wz"));
    if (bgzf_fp == nullptr) {
      throw std::runtime_error(fmt::format("Failed to open bgzip file for writing: {}", fname.string()));
    }
  } else {
    ofs.open(fname);
    if (!ofs) {
      throw std::runtime_error(fmt::format("Failed to open file for writing: {}", fname.string()));
    }
  }

  // Helper lambda to write to appropriate output
  auto write_line = [use_bgzip, &bgzf_fp, &fname, &ofs](const std::string& line) {
    const std::string data = line + "\n";
    if (use_bgzip) {
      if (bgzf_write(bgzf_fp.get(), data.c_str(), data.size()) < 0) {
        throw std::runtime_error(fmt::format("Failed to write to bgzip file: {}", fname.string()));
      }
    } else {
      ofs << line << '\n';
    }
  };

  // write the header to file
  write_line("##fileformat=VCFv4.2");
  write_line(fmt::format("##source={}", PROGRAM_NAME));
  if (!params.command_line_info.version.empty()) {
    write_line(io::GetCommandInfo(params.command_line_info));
  }
  write_line(fmt::format("##reference={}", params.reference_file));
  write_line(fmt::format("##Sex={}", enum_util::FormatEnumName(params.sex)));
  if (IsSomaticLikeMode(params.mode)) {
    WriteSomaticSampleHeaders(params, write_line);
  }
  if (!params.seq_lengths.empty()) {
    WriteContigHeaders(params, write_line);
  }

  WriteVcfMetadataLines(params, write_line);

  write_line(GetVcfHeaderLine(params.sample_id));
  for (const auto& seg : segments) {
    write_line(SegmentToRecordString(seg, params.mode));
  }

  // Close appropriate file handle
  if (use_bgzip) {
    bgzf_fp.reset();

    // Create tabix index for bgzipped VCF
    const s32 index_build_ret = tbx_index_build(fname.c_str(), 0, &tbx_conf_vcf);
    if (index_build_ret != 0) {
      throw std::runtime_error(fmt::format("Failed to create tabix index for: {}", fname.string()));
    }
  } else {
    ofs.close();
  }
}
}  // namespace xoos::cnc
