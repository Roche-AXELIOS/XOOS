#include "likelihood/flag-copy-number-calls.h"

#include <algorithm>
#include <cstddef>
#include <string>

#include <xoos/sex_predict/sex.h>
#include <xoos/stats/copy-number-stats.h>
#include <xoos/types/float.h>
#include <xoos/types/vec.h>

#include "likelihood/likelihood-flags.h"
#include "likelihood/likelihood-mode.h"
#include "likelihood/likelihood.h"
#include "segmentation/genomic-segments.h"

namespace xoos::cnc {

void FlagCallsByMeanMapq(vec<GenomicSegment>& segments, f64 mapq_cutoff_for_calls) {
  for (auto& segment : segments) {
    // skipped segments or segments with no MAPQ observations
    if (!segment.avg_mean_mapq.has_value()) {
      continue;
    }
    if ((!segment.total_copy_number.has_value() ||
         (segment.total_copy_number.has_value() && segment.total_copy_number.value() != 0)) &&
        segment.avg_mean_mapq.value() < mapq_cutoff_for_calls) {
      if (!segment.flags.has_value()) {
        segment.flags.emplace(vec<std::string>{});
      }
      segment.flags.value().emplace_back(kLikelihoodFlagLowAvgMeanMapq);
    }
  }
}

/**
 * @brief assigns the expected total copy number and then adds a flag if expected total copy number doesn't meet the
 * range for deletions or duplications. This should only be used in germline mode.
 * @param segments
 * @param sex
 * @param options
 */
void FlagCallsByExpectedTotalCopyNumber(vec<GenomicSegment>& segments,
                                        const sex_predict::Sex& sex,
                                        const NonIntegerTcnFilterOptions& options) {
  for (auto& seg : segments) {
    const bool expect_haploid = ExpectHaploid(seg, LikelihoodMode::kGermline, sex);
    const f64 mean_logr = seg.mean_logr.value();
    const f64 expected_total_copy_number =
        stats::ExpectedTotalCopyNumber(seg.purity.value(), seg.ploidy.value(), mean_logr, expect_haploid);
    seg.expected_total_copy_number = expected_total_copy_number;
    const bool het_del_not_in_range =
        *seg.total_copy_number == 1 && (expected_total_copy_number < options.het_del_min_expected_tcn ||
                                        expected_total_copy_number > options.het_del_max_expected_tcn);
    const bool hom_del_not_in_range =
        *seg.total_copy_number == 0 && expected_total_copy_number > options.hom_del_max_expected_tcn;
    const size_t reference_copy_number = expect_haploid ? 1 : kDiploidCopyNumber;
    const bool is_duplication = *seg.total_copy_number > reference_copy_number;
    bool dup_not_in_range = false;
    if (const auto& dup_min_expected_tcn = options.dup_min_expected_tcn;
        dup_min_expected_tcn.has_value() && is_duplication) {
      const f64 reference_adjusted_dup_min_expected_tcn =
          *dup_min_expected_tcn - static_cast<f64>(kDiploidCopyNumber) + static_cast<f64>(reference_copy_number);
      dup_not_in_range = expected_total_copy_number < reference_adjusted_dup_min_expected_tcn;
    }
    const auto& flag_to_assign = kLikelihoodFlagNonIntegerTcn;
    // remove the flag if it's already there
    if (seg.flags.has_value()) {
      auto& flags = seg.flags.value();
      if (std::ranges::find(flags, flag_to_assign) != flags.end()) {
        std::erase(flags, flag_to_assign);
      }
    }
    if (het_del_not_in_range || hom_del_not_in_range || dup_not_in_range) {
      if (!seg.flags.has_value()) {
        seg.flags.emplace(vec<std::string>{});
      }
      seg.flags.value().emplace_back(flag_to_assign);
    }
  }
}

void FlagCallsByLength(vec<GenomicSegment>& segments, size_t cnv_length_flag_min_size) {
  for (auto& seg : segments) {
    if (seg.end - seg.start < cnv_length_flag_min_size) {
      if (!seg.flags.has_value()) {
        seg.flags.emplace(vec<std::string>{});
      }
      seg.flags.value().emplace_back(kLikelihoodFlagCnvLength);
    }
  }
}

/**
 * @brief adds the kLikelihoodFlagChrYFemale flag to the segment "flag" members that are on chrY if the sample is female
 * @param segments
 */
void FlagChrYCallsIfFemale(vec<GenomicSegment>& segments) {
  for (auto& seg : segments) {
    if (sex_predict::IsInChromY(seg.contig) && seg.sex == sex_predict::Sex::kFemale) {
      if (!seg.flags.has_value()) {
        seg.flags.emplace(vec<std::string>{});
      }
      seg.flags.value().emplace_back(kLikelihoodFlagChrYFemale);
    }
  }
}

}  // namespace xoos::cnc
