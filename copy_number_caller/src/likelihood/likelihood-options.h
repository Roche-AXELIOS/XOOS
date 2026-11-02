#pragma once

#include <optional>

#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "likelihood/likelihood-mode.h"
#include "likelihood/likelihood-model.h"
#include "segmentation/segment-type.h"

namespace xoos::cnc {
const s32 kLikelihoodDefaultSomaticMAPQCutoffForCalls = -1;
const s32 kLikelihoodDefaultGermlineMAPQCutoffForCalls = 30;
const size_t kLikelihoodCnvLengthFlagMinSize = 2000;
const f64 kLikelihoodDefaultHetDelMinExpectedTcn = 0.7;
const f64 kLikelihoodDefaultHetDelMaxExpectedTcn = 1.3;
const f64 kLikelihoodDefaultHomDelMaxExpectedTcn = 0.3;
const f64 kLikelihoodSbxSHetDelMinExpectedTcn = 0.8;
const f64 kLikelihoodSbxSHetDelMaxExpectedTcn = 1.2;
const f64 kLikelihoodSbxSDupMinExpectedTcn = 2.65;

struct NonIntegerTcnFilterOptions {
  f64 het_del_min_expected_tcn{kLikelihoodDefaultHetDelMinExpectedTcn};
  f64 het_del_max_expected_tcn{kLikelihoodDefaultHetDelMaxExpectedTcn};
  f64 hom_del_max_expected_tcn{kLikelihoodDefaultHomDelMaxExpectedTcn};
  // Only SBX-S germline has a calibrated duplication threshold; other library types leave this unset.
  std::optional<f64> dup_min_expected_tcn;
};

const NonIntegerTcnFilterOptions kLikelihoodDefaultNonIntegerTcnFilterOptions{};
const NonIntegerTcnFilterOptions kLikelihoodSbxSNonIntegerTcnFilterOptions{
    .het_del_min_expected_tcn = kLikelihoodSbxSHetDelMinExpectedTcn,
    .het_del_max_expected_tcn = kLikelihoodSbxSHetDelMaxExpectedTcn,
    .hom_del_max_expected_tcn = kLikelihoodDefaultHomDelMaxExpectedTcn,
    .dup_min_expected_tcn = kLikelihoodSbxSDupMinExpectedTcn,
};

struct LikelihoodOptions {
  LikelihoodModel likelihood_model = LikelihoodModel::kSerialSummarized;
  SegmentType input_segments_type{SegmentType::kUnknown};
  SegmentType output_segments_type{SegmentType::kUnknown};
  size_t cnv_length_flag_min_size = kLikelihoodCnvLengthFlagMinSize;
  s32 mapq_cutoff_for_calls{};
  bool mapq_cutoff_for_calls_is_user_set{false};
  NonIntegerTcnFilterOptions non_integer_tcn_filter_options{};
  LikelihoodMode mode{LikelihoodMode::kUnknown};
  // Tumor-only mode: after somatic prediction, classify each segment as germline or somatic and emit the
  // IsGermlineCnv SEG column and VCF INFO/GERMLINE flag.
  bool tumor_only{false};
};

}  // namespace xoos::cnc
