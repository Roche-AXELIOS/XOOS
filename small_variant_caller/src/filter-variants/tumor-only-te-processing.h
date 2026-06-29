#pragma once

#include <optional>
#include <string>

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/filtering.h"
#include "core/score-calculator.h"
#include "core/variant-id.h"
#include "core/variant-info.h"
#include "util/parallel-compute-utils.h"

namespace xoos::svc {

class TumorOnlyTeProcessing {
 public:
  TumorOnlyTeProcessing(const GlobalContext& global_ctx,
                        const FilterSettings& filter_settings,
                        const PhasedFilterSettings& phased_filter_settings,
                        u64 prev_pos,
                        VarIdToVcfFeatures vcf_features,
                        BamRegionFeatureCollection bam_features);

  /**
   * @brief Score a variant using the somatic model and update its filter status. The variant score is set directly in
   * the internal BAM feature map for use later.
   * @param vid The VariantId of the variant to score
   * @param normalize_targets A map of normalize targets for each chromosomes
   * @param scoring_cols The feature columns to use for scoring
   * @param somatic_calculator The ScoreCalculator instance to use for scoring
   * @return An optional BamFeatureTuple if the variant was found and scored, otherwise std::nullopt
   */
  std::optional<BamFeatureTuple> ScoreVariant(const VariantId& vid,
                                              const ChromMedianDepth& normalize_targets,
                                              const vec<FeatureColumn>& scoring_cols,
                                              const ScoreCalculator& somatic_calculator);

  /**
   * @brief Create new VCF records for any forced calls between the previous position and the current position.
   * @param pos The current position being processed
   * @param chrom The chromosome the current position is on
   * @param out_records Output vector to store newly created VCF records
   */
  void CreateNewForceRecordsForInbetweenPositions(u64 pos,
                                                  const std::string& chrom,
                                                  vec<io::VcfRecordPtr>& out_records);

  /**
   * @brief Process a VCF record, updating it with features and filter status as needed. If the record is phased, a
   * different set of filtering criteria is applied.
   * @param vid The VariantId of the record
   * @param record A VCF record to process
   * @param bam_feature An optional UnifiedVariantFeature if available for the variant
   * @param out_records A vector to store processed VCF records for output
   * @param make_forcedcalls True if forced calls need to be made
   * @param key A variant correlation key for the variant
   * @param key_unpadded An unpadded variant correlation key for the variant
   */
  void ProcessRecord(const VariantId& vid,
                     const io::VcfRecordPtr& record,
                     const std::optional<BamFeatureTuple>& bam_feature,
                     vec<io::VcfRecordPtr>& out_records,
                     bool make_forcedcalls,
                     const std::string& key,
                     const std::string& key_unpadded);

 private:
  const FilterSettings& _filter_settings;
  const PhasedFilterSettings& _phased_filter_settings;
  u64 _prev_pos;
  VarIdToVcfFeatures _vcf_features;
  BamRegionFeatureCollection _bam_features;
  io::VcfHeaderPtr _hdr;
  StrUnorderedMap<std::unordered_set<u64>> _forced_call_chrom_pos;
  StrUnorderedMap<std::unordered_map<u64, vec<VariantId>>> _forced_call_chrom_pos_to_vids;
  bool _phased;

  StrUnorderedMap<std::unordered_set<u64>> _seen_forced_call_chrom_pos{};
};

}  // namespace xoos::svc
