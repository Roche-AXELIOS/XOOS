#pragma once

#include <optional>

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "compute-vcf-features/compute-vcf-features.h"
#include "core/bam-feature-collection.h"
#include "core/score-calculator.h"
#include "core/variant-info.h"
#include "util/parallel-compute-utils.h"

namespace xoos::svc {

const s32 kTumorNormalGq{100};

class TumorNormalProcessing {
 public:
  TumorNormalProcessing(const GlobalContext& global_ctx,
                        VarIdToVcfFeatures vcf_features,
                        BamRegionFeatureCollection bam_features);

  /**
   * @brief Process a Tumor Normal VCF record, updating it with features and filter status as needed.
   * @param record A VCF record to process
   * @param out_records A vector to store processed VCF records for output
   * @param somatic_calculator A ScoreCalculator object for somatic variant scoring
   * @param normalize_target A DpTuple for feature normalization
   * @param header_info A VcfHeaderInfo object containing VCF header metadata
   * @param shap_value_rows A vector to store SHAP value rows for output, if applicable
   */
  void ProcessRecord(const io::VcfRecordPtr& record,
                     vec<io::VcfRecordPtr>& out_records,
                     const ScoreCalculator& somatic_calculator,
                     const DepthTuple& normalize_target,
                     const VcfHeaderInfo& header_info,
                     vec<vec<std::string>>& shap_value_rows) const;

 private:
  std::string_view CheckHardFilters(const TumorNormalBamFeatureTuple& bam_feat,
                                    const VcfFeature& vcf_feat,
                                    const DepthTuple& normalize_target) const;

  std::pair<std::string_view, PredictionScore> ScoreVariantML(const VariantId& vid,
                                                              const TumorNormalBamFeatureTuple& bam_feat,
                                                              const VcfFeature& vcf_feat,
                                                              const DepthTuple& normalize_target,
                                                              const ScoreCalculator& calculator,
                                                              vec<vec<std::string>>& shap_value_rows) const;

  VarIdToVcfFeatures _vcf_features;
  BamRegionFeatureCollection _bam_features;
  io::VcfHeaderPtr _hdr;
  f32 _snv_min_ml_score;
  f32 _indel_min_ml_score;
  u32 _min_tumor_support;
  u32 _max_normal_support;
  f32 _min_tumor_af;
  // Relaxed thresholds applied when a variant is flagged `in_catalog` (--catalog-vcf).
  // When unset, the corresponding standard threshold above is used (no relaxation).
  std::optional<f32> _catalog_min_tumor_af;
  std::optional<f32> _catalog_snv_min_ml_score;
  std::optional<f32> _catalog_indel_min_ml_score;
  f32 _min_dp_ratio;
  u32 _max_indel_size;
  bool _normalize_scoring_features;
  vec<FeatureColumn> _scoring_cols;
  bool _is_duplex_protocol;
};
}  // namespace xoos::svc
