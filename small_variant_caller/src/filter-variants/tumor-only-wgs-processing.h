#pragma once

#include <utility>

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/bam-feature-collection.h"
#include "core/score-calculator.h"
#include "core/variant-info.h"
#include "util/parallel-compute-utils.h"

namespace xoos::svc {

/**
 * @brief Per-region processing for the `tumor-only-wgs` workflow.
 *
 * Scores single-sample variants with a single binary somatic model and applies pre-ML hard filters. Unlike
 * `tumor-normal-wgs`, there is no normal sample: features are read from the unprefixed single-sample BAM feature
 * tuple, and the only hard filter is the minimum alt-support gate. SNV and indel variants use separate ML score
 * thresholds against the same model. Panel of Normals features are consumed transparently via the scoring columns
 * and default to 0 when no PON VCF was supplied to `compute_vcf_features`.
 */
class TumorOnlyWgsProcessing {
 public:
  TumorOnlyWgsProcessing(const GlobalContext& global_ctx,
                         VarIdToVcfFeatures vcf_features,
                         BamRegionFeatureCollection bam_features);

  /**
   * @brief Process a single tumor-only-wgs VCF record, updating it with features and filter status.
   * @param record A VCF record to process
   * @param out_records A vector to store processed VCF records for output
   * @param somatic_calculator A ScoreCalculator object for somatic variant scoring
   * @param normalize_target A DepthTuple for feature normalization
   * @param shap_value_rows A vector to store SHAP value rows for output
   */
  void ProcessRecord(const io::VcfRecordPtr& record,
                     vec<io::VcfRecordPtr>& out_records,
                     const ScoreCalculator& somatic_calculator,
                     const DepthTuple& normalize_target,
                     vec<vec<std::string>>& shap_value_rows) const;

 private:
  /**
   * @brief Check pre-ML hard filters for a single-sample variant.
   * @param bam_feat BAM feature tuple for the variant.
   * @return Filter ID of the first failing threshold, or `kFilteringPassId`.
   */
  std::string_view CheckHardFilters(const BamFeatureTuple& bam_feat) const;

  /**
   * @brief Run the somatic ML model, record SHAP values, and return the resulting filter and score.
   * @param vid Variant identifier.
   * @param bam_feat BAM feature tuple for the variant.
   * @param vcf_feat VCF feature values for the variant.
   * @param normalize_target Chromosome median depth targets.
   * @param calculator Score calculator wrapping the ML model.
   * @param shap_value_rows Output accumulator for SHAP value rows.
   * @return Pair of {filter ID, prediction score}.
   */
  std::pair<std::string_view, PredictionScore> ScoreVariantML(const VariantId& vid,
                                                              const BamFeatureTuple& bam_feat,
                                                              const VcfFeature& vcf_feat,
                                                              const DepthTuple& normalize_target,
                                                              const ScoreCalculator& calculator,
                                                              vec<vec<std::string>>& shap_value_rows) const;

  VarIdToVcfFeatures _vcf_features;
  BamRegionFeatureCollection _bam_features;
  io::VcfHeaderPtr _hdr;
  f32 _snv_min_ml_score;
  f32 _indel_min_ml_score;
  u32 _min_support;
  f32 _min_af;
  u32 _max_indel_size;
  bool _normalize_scoring_features;
  vec<FeatureColumn> _scoring_cols;
  bool _is_duplex_protocol;
};

}  // namespace xoos::svc
