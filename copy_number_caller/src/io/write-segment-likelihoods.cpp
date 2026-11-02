#include "io/write-segment-likelihoods.h"

#include <concepts>
#include <stdexcept>

#include <xoos/enum/enum-util.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/util/missing-value.h>

#include "io/column-names.h"
#include "util/parse-util.h"

namespace xoos::cnc {

namespace {

constexpr s32 kPurityPloidyOutputPrecision = 2;

// Streams an optional value to the TSV, emitting the standard missing-value token when the value is absent. Floating
// point values go through FloatAsString for fixed 6-decimal-place output; other types stream directly.
template <typename T>
void WriteOptional(std::ofstream& ofs, const std::optional<T>& value) {
  if (value.has_value()) {
    if constexpr (std::floating_point<T>) {
      ofs << FloatAsString(value.value());
    } else {
      ofs << value.value();
    }
  } else {
    ofs << util::kNA;
  }
}

}  // namespace

SegmentLikelihoodWriter::SegmentLikelihoodWriter(const SegmentLikelihoodsParams& params)
    : _ofs{params.out_fname}, _tumor_only{params.tumor_only} {
  if (!_ofs.is_open()) {
    Logging::Error("Failed to open segment likelihood output file: {}", params.out_fname.string());
    throw std::runtime_error("Failed to open segment likelihood output file: " + params.out_fname.string());
  }
  if (!params.command_line_info.version.empty()) {
    io::WriteTsvMetadata(_ofs, params.command_line_info);
  }
  // Record which purity/ploidy search pass produced this file so the two per-pass diagnostics are self-describing.
  _ofs << "##TotalCopyNumberPrior=" << enum_util::FormatEnumName(params.prior) << "\n";

  // Fixed leading columns, then the somatic-state MLE columns. The tumor-only layout appends germline likelihood/TCN,
  // per-cell classification, and contribution likelihood; the tumor-normal layout omits them.
  _ofs << kColumnContig << "\t" << kColumnStart << "\t" << kColumnEnd << "\t" << kSegmentLikelihoodsOutPurity << "\t"
       << kSegmentLikelihoodsOutPloidy << "\t" << kSegmentLikelihoodsOutSomaticLogRLikelihood << "\t"
       << kSegmentLikelihoodsOutSomaticBafLikelihood << "\t" << kSegmentLikelihoodsOutSomaticJointLikelihood << "\t"
       << kSegmentLikelihoodsOutSomaticTotalCopyNumber << "\t" << kSegmentLikelihoodsOutSomaticMinorCopyNumber;
  if (_tumor_only) {
    _ofs << "\t" << kSegmentLikelihoodsOutGermlineLogRLikelihood << "\t"
         << kSegmentLikelihoodsOutGermlineTotalCopyNumber << "\t" << kSegmentLikelihoodsOutClassification << "\t"
         << kSegmentLikelihoodsOutContributionLikelihood;
  }
  _ofs << "\n";
}

void SegmentLikelihoodWriter::WriteRow(const SegmentLikelihoodRowKey& key,
                                       const SomaticLikelihoodRow& somatic,
                                       const std::optional<GermlineLikelihoodRow>& germline,
                                       const std::optional<bool> classified_germline) {
  // Purity/ploidy match the two-decimal-place format used by related outputs, while likelihoods use six decimal places.
  _ofs << key.contig << "\t" << key.start << "\t" << key.end << "\t"
       << FloatAsString(key.purity, kPurityPloidyOutputPrecision) << "\t"
       << FloatAsString(key.ploidy, kPurityPloidyOutputPrecision) << "\t" << FloatAsString(somatic.logr_likelihood)
       << "\t";
  WriteOptional(_ofs, somatic.baf_likelihood);
  _ofs << "\t";
  WriteOptional(_ofs, somatic.joint_likelihood);
  _ofs << "\t" << somatic.total_copy_number << "\t";
  WriteOptional(_ofs, somatic.minor_copy_number);

  if (_tumor_only) {
    _ofs << "\t";
    if (germline.has_value()) {
      _ofs << FloatAsString(germline.value().logr_likelihood) << "\t" << germline.value().total_copy_number;
    } else {
      _ofs << util::kNA << "\t" << util::kNA;
    }
    _ofs << "\t";
    if (classified_germline.has_value()) {
      const bool is_germline = classified_germline.value();
      // The contribution likelihood is what the segment actually contributes to the grid sum at this cell: its germline
      // logR likelihood when classified germline, otherwise its somatic joint likelihood.
      const f64 contribution =
          is_germline ? germline.value().logr_likelihood : somatic.joint_likelihood.value_or(somatic.logr_likelihood);
      _ofs << (is_germline ? kSegmentClassificationGermline : kSegmentClassificationSomatic) << "\t"
           << FloatAsString(contribution);
    } else {
      _ofs << util::kNA << "\t" << util::kNA;
    }
  }
  _ofs << "\n";
}

}  // namespace xoos::cnc
