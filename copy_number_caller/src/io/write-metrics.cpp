#include "io/write-metrics.h"

#include <cstddef>
#include <fstream>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>

#include <fmt/core.h>
#include <fmt/format.h>

#include <xoos/enum/enum-util.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/vec.h>
#include <xoos/util/missing-value.h>

#include "cli/copy-number-caller-modes.h"
#include "coverage/coverage-check.h"
#include "intervals/augmented-intervals.h"
#include "likelihood/likelihood-flags.h"
#include "observations/vcf-check.h"
#include "seg-to-vcf/seg-to-vcf.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {

namespace {

/// Per-status segment tallies, each paired with the subset that passed filtering.
struct SegmentCounts {
  std::size_t gain{0};
  std::size_t loss{0};
  std::size_t cnloh{0};
  std::size_t gainloh{0};
  std::size_t passing_gain{0};
  std::size_t passing_loss{0};
  std::size_t passing_cnloh{0};
  std::size_t passing_gainloh{0};
};

/// Everything the row writers need, so each takes a stream and a single summary.
struct MetricsSummary {
  std::size_t num_targets{0};
  std::size_t n_snps{0};
  sex_predict::Sex sex{};
  std::size_t num_segments{0};
  std::optional<f64> purity;
  std::optional<f64> ploidy;
  SegmentCounts counts;
};

/// Germline and somatic modes classify segments differently; other modes contribute no counts.
CopyNumberStatus ClassifySegment(const GenomicSegment& seg, const CopyNumberCallerModes mode) {
  if (mode == CopyNumberCallerModes::kGermlineNormalWGS) {
    return GetCopyNumberStatus(seg);
  }
  if (IsSomaticLikeMode(mode)) {
    return GetSomaticCopyNumberStatus(seg);
  }
  return CopyNumberStatus::kUnknown;
}

SegmentCounts CountSegmentsByStatus(const vec<GenomicSegment>& segments, const CopyNumberCallerModes mode) {
  using enum CopyNumberStatus;
  SegmentCounts counts;
  for (const auto& seg : segments) {
    const std::size_t passing =
        ConstructVcfFilterFieldFromGenomicSegment(seg) == kLikelihoodFlagPass ? std::size_t{1} : std::size_t{0};
    switch (ClassifySegment(seg, mode)) {
      case kGain:
        ++counts.gain;
        counts.passing_gain += passing;
        break;
      case kLoss:
        ++counts.loss;
        counts.passing_loss += passing;
        break;
      case kCnLoh:
        ++counts.cnloh;
        counts.passing_cnloh += passing;
        break;
      case kGainLoh:
        ++counts.gainloh;
        counts.passing_gainloh += passing;
        break;
      default:
        break;
    }
  }
  return counts;
}

void WriteGermlineRows(std::ostream& ofs, const MetricsSummary& summary) {
  ofs << "Metric\tValue\n"
      << "Number of targets\t" << summary.num_targets << '\n'
      << "Sex\t" << enum_util::FormatEnumName(summary.sex) << '\n'
      << "Number of segments\t" << summary.num_segments << '\n'
      << "Number of gain segments\t" << summary.counts.gain << '\n'
      << "Number of loss segments\t" << summary.counts.loss << '\n'
      << "Number of passing gain segments\t" << summary.counts.passing_gain << '\n'
      << "Number of passing loss segments\t" << summary.counts.passing_loss << '\n';
}

void WriteSomaticRows(std::ostream& ofs, const MetricsSummary& summary) {
  const std::string purity_str = util::FormatOrNA(summary.purity, "{:.2f}");
  const std::string ploidy_str = util::FormatOrNA(summary.ploidy, "{:.2f}");
  ofs << "Metric\tValue\n"
      << "Number of targets\t" << summary.num_targets << '\n'
      << "Number of SNPs\t" << summary.n_snps << '\n'
      << "Sex\t" << enum_util::FormatEnumName(summary.sex) << '\n'
      << "Tumor Purity\t" << purity_str << '\n'
      << "Tumor Ploidy\t" << ploidy_str << '\n'
      << "Number of segments\t" << summary.num_segments << '\n'
      << "Number of somatic gain (non-LOH) segments\t" << summary.counts.gain << '\n'
      << "Number of somatic loss segments\t" << summary.counts.loss << '\n'
      << "Number of somatic copy-neutral LOH segments\t" << summary.counts.cnloh << '\n'
      << "Number of somatic gain LOH segments\t" << summary.counts.gainloh << '\n'
      << "Number of passing somatic gain (non-LOH) segments\t" << summary.counts.passing_gain << '\n'
      << "Number of passing somatic loss segments\t" << summary.counts.passing_loss << '\n'
      << "Number of passing somatic copy-neutral LOH segments\t" << summary.counts.passing_cnloh << '\n'
      << "Number of passing somatic gain LOH segments\t" << summary.counts.passing_gainloh << '\n';
}

void WriteCoverageRows(std::ostream& ofs, const CoverageCheckResult& coverage_check) {
  ofs << fmt::format("Median coverage\t{:.2f}\n", coverage_check.median_coverage)
      << fmt::format("Percentage of low coverage windows\t{:.2f}\n", coverage_check.pct_windows_below_threshold);
}

/// Input-VCF metrics that only a matched normal can produce, so they are tumor-normal only.
void WriteTumorNormalVcfRows(std::ostream& ofs, const VcfCheckResult& vcf_check) {
  const std::string het_snp_fraction_str = util::FormatOrNA(vcf_check.het_snp_fraction, "{:.2f}");
  ofs << "Number of heterozygous germline SNPs\t" << vcf_check.num_het_snps << '\n'
      << "Number of somatic variants\t" << vcf_check.num_somatic_variants << '\n'
      << "Het SNP fraction\t" << het_snp_fraction_str << '\n';
}

}  // namespace

void WriteMetricsFile(const fs::path& metrics_fname,
                      const AugmentedIntervals& augmented_intervals,
                      const size_t n_snps,
                      const sex_predict::Sex sex,
                      const vec<GenomicSegment>& segments,
                      const std::optional<f64>& purity,
                      const std::optional<f64>& ploidy,
                      const CopyNumberCallerModes mode,
                      const io::CommandLineInfo& command_line_info,
                      const std::optional<CoverageCheckResult>& coverage_check,
                      const std::optional<VcfCheckResult>& vcf_check) {
  std::ofstream metrics_ofs(metrics_fname);
  if (!metrics_ofs.is_open()) {
    Logging::Error("Failed to open metrics file: {}", metrics_fname.string());
    throw std::runtime_error("Failed to open metrics file: " + metrics_fname.string());
  }
  if (!command_line_info.version.empty()) {
    io::WriteTsvMetadata(metrics_ofs, command_line_info);
  }

  const MetricsSummary summary{.num_targets = augmented_intervals.GetRegions().size(),
                               .n_snps = n_snps,
                               .sex = sex,
                               .num_segments = segments.size(),
                               .purity = purity,
                               .ploidy = ploidy,
                               .counts = CountSegmentsByStatus(segments, mode)};

  if (mode == CopyNumberCallerModes::kGermlineNormalWGS) {
    WriteGermlineRows(metrics_ofs, summary);
  } else if (IsSomaticLikeMode(mode)) {
    WriteSomaticRows(metrics_ofs, summary);
  }

  if (coverage_check.has_value()) {
    WriteCoverageRows(metrics_ofs, *coverage_check);
  }

  // Gated on the mode, not just on the data: tumor-only forwards its (always empty) vcf_check into
  // this writer, so relying on the caller alone would let these rows leak into tumor-only output the
  // moment that field is ever populated.
  if (mode == CopyNumberCallerModes::kSomaticTumorNormalWGS && vcf_check.has_value()) {
    WriteTumorNormalVcfRows(metrics_ofs, *vcf_check);
  }
}

}  // namespace xoos::cnc
