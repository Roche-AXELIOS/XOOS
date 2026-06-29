#pragma once

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
                        u64 prev_pos,
                        VarIdToVcfFeatures vcf_features,
                        BamRegionFeatureCollection bam_features);

  /**
   * @brief Process a Tumor Normal VCF record, updating it with features and filter status as needed.
   * @param record A VCF record to process
   * @param out_records A vector to store processed VCF records for output
   * @param somatic_calculator A ScoreCalculator object for somatic variant scoring
   * @param normalize_target A DpTuple for feature normalization
   * @param header_info A VcfHeaderInfo object containing VCF header metadata
   */
  void ProcessRecord(const io::VcfRecordPtr& record,
                     vec<io::VcfRecordPtr>& out_records,
                     const ScoreCalculator& somatic_calculator,
                     const DepthTuple& normalize_target,
                     const VcfHeaderInfo& header_info) const;

 private:
  u64 _prev_pos;
  VarIdToVcfFeatures _vcf_features;
  BamRegionFeatureCollection _bam_features;
  io::VcfHeaderPtr _hdr;
  f32 _somatic_ml_snv_threshold;
  f32 _somatic_ml_indel_threshold;
  u32 _tumor_support_threshold;
  vec<FeatureColumn> _scoring_cols;
};
}  // namespace xoos::svc
