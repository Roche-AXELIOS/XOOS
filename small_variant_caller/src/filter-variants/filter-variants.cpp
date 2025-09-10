#include "filter-variants.h"

#include <algorithm>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>

#include "compute-bam-features/progress-meter.h"
#ifdef SOMATIC_ENABLE
#include <fstream>
#include <sstream>
#include <unordered_set>
#endif  // SOMATIC_ENABLE

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/log/logging.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/string-functions.h>

#include "compute-vcf-features/compute-vcf-features.h"
#include "core/command-line-info.h"
#include "core/filtering.h"
#include "core/model-metadata.h"
#include "core/vcf-fields.h"
#include "util/log-util.h"
#include "util/parallel-compute-utils.h"

namespace xoos::svc {

using BedRegion = io::BedRegion;

#ifdef SOMATIC_ENABLE
/**
 * Loads a block list from file
 * @param block_list The blocklist in chr_pos_ref_alt format
 * @return The parsed left-padded blocklist
 */
StrUnorderedSet LoadBlockList(const fs::path& block_list) {
  StrUnorderedSet result;
  std::ifstream block_list_stream(block_list);
  std::string line;
  while (getline(block_list_stream, line)) {
    // Need to left pad the blocklist so it matches the behavior of GetVariantCorrelationKey
    vec<std::string> split_line;
    std::string seg;
    std::stringstream input_line(line);
    while (std::getline(input_line, seg, '_')) {
      split_line.push_back(seg);
    }
    if (split_line.size() == 4) {
      std::string key =
          GetVariantCorrelationKey(split_line[0], std::stoul(split_line[1]) - 1, split_line[2], split_line[3]);
      result.insert(key);
    }
  }
  return result;
}

/**
 * @brief Update a VCF record with selected feature information. Intended for `somatic` workflow only.
 * @param record VCF record to be updated
 * @param var_feat Variant feature struct
 * @param ref_feat Reference feature struct
 * @param is_hotspot Indicates if the variant is a hotspot
 */
void UpdateRecord(const io::VcfRecordPtr& record,
                  const UnifiedVariantFeature& var_feat,
                  const UnifiedReferenceFeature& ref_feat,
                  const bool is_hotspot) {
  record->SetFormatField(kWeightedCountsId, vec<f32>{static_cast<f32>(var_feat.weighted_score)});
  record->SetFormatField(kMachineLearningId, vec<f32>{static_cast<f32>(var_feat.ml_score)});
  record->SetFormatField(kNonDuplexCountsId, vec<s32>{static_cast<s32>(var_feat.nonduplex)});
  record->SetFormatField(kDuplexCountsId, vec<s32>{static_cast<s32>(var_feat.duplex)});
  record->SetFormatField(kStrandBiasMetricId, vec<f32>{static_cast<f32>(var_feat.strandbias)});
  record->SetFormatField(kPlusOnlyCountsId, vec<s32>{static_cast<s32>(var_feat.plusonly)});
  record->SetFormatField(kMinusOnlyCountsId, vec<s32>{static_cast<s32>(var_feat.minusonly)});
  record->SetFormatField(kSequenceContextId, vec<std::string>{var_feat.context});
  // update the existing AD field, since we use different read filtering than Mutect2
  record->SetFormatField(
      kAlleleDepthId,
      vec<s32>{static_cast<s32>(ref_feat.support), static_cast<s32>(var_feat.duplex + var_feat.nonduplex)});
  if (is_hotspot) {
    record->SetFormatField(kHotspotId, vec<s32>{1});
  }
}

/**
 * @brief Update a phased VCF record with empty values. Intended for `somatic` workflow only.
 * @param record VCF record to be updated
 */
void UpdatePhasedRecord(const io::VcfRecordPtr& record) {
  // Phased calls don't use ML based filtering and do not have values for these fields
  record->SetFormatField(kWeightedCountsId, vec<f32>{static_cast<f32>(0.0)});
  record->SetFormatField(kMachineLearningId, vec<f32>{static_cast<f32>(0.0)});
  record->SetFormatField(kNonDuplexCountsId, vec<s32>{0});
  record->SetFormatField(kDuplexCountsId, vec<s32>{0});
  record->SetFormatField(kStrandBiasMetricId, vec<f32>{static_cast<f32>(0.0)});
  record->SetFormatField(kPlusOnlyCountsId, vec<s32>{0});
  record->SetFormatField(kMinusOnlyCountsId, vec<s32>{0});
  record->SetFormatField(kSequenceContextId, vec<std::string>{"NA"});
}

/**
 * @brief Set up a VCF record for forced variant calls.
 * @param record VCF record to be updated
 * @param variant_id VariantId struct
 */
void SetupForceRecord(const io::VcfRecordPtr& record, const VariantId& variant_id) {
  const vec<std::string> alleles = {variant_id.ref, variant_id.alt};
  record->SetAlleles(alleles);
  record->SetFilter(kFilteringForcedId);
}

void SetSomaticTNValues(const io::VcfRecordPtr& record,
                        const VariantId& variant_id,
                        const UnifiedVariantFeature& bam_feat,
                        const UnifiedReferenceFeature& ref_feat,
                        const VcfFeature& vcf_feat,
                        const bool tumor_first,
                        const f64 pred_score) {
  // Setup the format and info fields for tumor normal vcf records. Assumes only two samples in the VCF, one tumor, one
  // normal
  const auto normal_gt = "0/0";
  const auto tumor_gt = "0/1";
  record->SetFormatField(kFieldGt,
                         tumor_first ? vec<std::string>{tumor_gt, normal_gt} : vec<std::string>{normal_gt, tumor_gt});
  const auto tumor_support = static_cast<s32>(bam_feat.tumor_support);
  const auto normal_support = static_cast<s32>(bam_feat.normal_support);
  const auto tumor_ref_support = static_cast<s32>(ref_feat.tumor_support);
  const auto normal_ref_support = static_cast<s32>(ref_feat.normal_support);
  record->SetFormatField(kFieldAd,
                         tumor_first ? vec<s32>{tumor_ref_support, tumor_support, normal_ref_support, normal_support}
                                     : vec<s32>{normal_ref_support, normal_support, tumor_ref_support, tumor_support});
  const auto tumor_af = static_cast<f32>(bam_feat.tumor_af);
  const auto normal_af = static_cast<f32>(bam_feat.normal_af);
  record->SetFormatField(kFieldAf, tumor_first ? vec<f32>{tumor_af, normal_af} : vec<f32>{normal_af, tumor_af});
  record->SetFormatField(kFieldDp,
                         tumor_first
                             ? vec<s32>{tumor_support + tumor_ref_support, normal_support + normal_ref_support}
                             : vec<s32>{normal_support + normal_ref_support, tumor_support + tumor_ref_support});
  const auto tumor_bq = static_cast<f32>(bam_feat.tumor_baseq_mean);
  const auto normal_bq = static_cast<f32>(bam_feat.normal_baseq_mean);
  record->SetFormatField(kBaseqQualId, tumor_first ? vec<f32>{tumor_bq, normal_bq} : vec<f32>{normal_bq, tumor_bq});
  const auto tumor_mq = static_cast<f32>(bam_feat.tumor_mapq_mean);
  const auto normal_mq = static_cast<f32>(bam_feat.normal_mapq_mean);
  record->SetFormatField(kMapQualId, tumor_first ? vec<f32>{tumor_mq, normal_mq} : vec<f32>{normal_mq, tumor_mq});
  const auto tumor_distance = static_cast<f32>(bam_feat.tumor_distance_mean);
  const auto normal_distance = static_cast<f32>(bam_feat.normal_distance_mean);
  record->SetFormatField(
      kDistanceId, tumor_first ? vec<f32>{tumor_distance, normal_distance} : vec<f32>{normal_distance, tumor_distance});
  record->SetFormatField(kFieldGq, vec<s32>{100, 100});
  record->SetInfoField<f32>(kPredId, {static_cast<f32>(pred_score)});
  record->SetInfoField<f32>(kRefBQId, {static_cast<f32>(ref_feat.baseq_mean)});
  record->SetInfoField<f32>(kRefMQId, {static_cast<f32>(ref_feat.mapq_mean)});
  record->SetInfoField<f32>(kAltBQId, {static_cast<f32>(bam_feat.baseq_mean)});
  record->SetInfoField<f32>(kAltMQId, {static_cast<f32>(bam_feat.mapq_mean)});
  record->SetInfoField<f32>(kFieldNalod, {vcf_feat.nalod});
  record->SetInfoField<f32>(kFieldNlod, {vcf_feat.nlod});
  record->SetInfoField<f32>(kFieldTlod, {vcf_feat.tlod});
  record->SetInfoField<s32>(kFieldMpos, {static_cast<s32>(vcf_feat.mpos)});
  record->SetInfoField<s32>(kSubtypeId, {variant_id.sub_index});
  record->SetInfoField<std::string>(kContextId, {bam_feat.context});
  record->SetInfoField<f32>(kFieldPopaf, {vcf_feat.popaf});
}

/**
 * @brief Copy a VCF record for a given allele within the somaticTN workflow.
 * @param original_record VCF record to be copied
 * @param new_header VCF header for the new VCF record
 * @return The copied VCF record
 */
io::VcfRecordPtr CopySomaticRecord(const io::VcfRecordPtr& original_record, const io::VcfHeaderPtr& new_header) {
  // TODO : Update this function to output required info
  const auto& record_copy = original_record->Clone(new_header);
  record_copy->SetAlleles({original_record->Allele(0), original_record->Allele(0)});
  return record_copy;
}

/**
 * @brief Update a VCF record for the given features within the somatic-tumor-normal workflow.
 * @param record VCF record to be updated
 * @param bam_feat BAM features
 * @param vcf_feat VCF features
 * @param ref_feat Reference features
 * @param score Genotype score
 */
void UpdateSomaticTNRecord(const io::VcfRecordPtr& record,
                           const UnifiedVariantFeature& bam_feat,
                           const VcfFeature& vcf_feat,
                           const UnifiedReferenceFeature& ref_feat,
                           const Genotype score) {
  UpdateGenotypeRelatedFields(record, score);
  record->SetFormatField<s32>(kGermlineMLId, {1});
  record->SetFormatField<s32>(kAlleleDepthId, {static_cast<s32>(ref_feat.support), static_cast<s32>(bam_feat.support)});
  record->SetFormatField<s32>(kFieldDp, {static_cast<s32>(ref_feat.support + bam_feat.support)});
  record->SetFormatField<f32>(kGermlineGnomadAFId, {vcf_feat.popaf});
  record->SetFormatField<f32>(kGermlineRefAvgMapqId, {static_cast<f32>(ref_feat.mapq_mean)});
  record->SetFormatField<f32>(kGermlineAltAvgMapqId, {static_cast<f32>(bam_feat.mapq_mean)});
  record->SetFormatField<f32>(kGermlineRefAvgDistId, {static_cast<f32>(ref_feat.distance_mean)});
  record->SetFormatField<f32>(kGermlineAltAvgDistId, {static_cast<f32>(bam_feat.distance_mean)});
  record->SetFormatField<f32>(kGermlineDensity100BPId, {static_cast<f32>(vcf_feat.variant_density)});
}

/**
 * @brief Fail a VCF record and update relevant fields.
 * @param record VCF record
 */
void FailSomaticTNRecord(const io::VcfRecordPtr& record, const std::string& fail_id) {
  record->AddFilter(fail_id);
}
#endif  // SOMATIC_ENABLE

/**
 * @brief Helper function to extract sex chromosome name and PARs from a BED file.
 * @param bed_path Path of BED file containing PARs
 * @return Chromosome name and vector of PARs.
 */
static std::pair<std::string, vec<Interval>> ExtractSexChromPAR(const std::optional<fs::path>& bed_path,
                                                                const std::string& chrom) {
  std::string name;
  vec<Interval> par;
  if (bed_path.has_value()) {
    auto chr_to_intervals = GetChromIntervalMap(bed_path).value();
    if (chr_to_intervals.empty()) {
      WarnAsErrorIfSet("No intervals found in {}", bed_path.value());
    } else {
      if (chr_to_intervals.size() == 1) {
        name = chr_to_intervals.begin()->first;
        par = chr_to_intervals.begin()->second;
      } else {
        vec<std::string> ref_names;
        for (const auto& [chrom, intervals] : chr_to_intervals) {
          ref_names.emplace_back(chrom);
        }
        WarnAsErrorIfSet("PARs BED file has multiple reference names: {}", string::Join(ref_names, ","));
      }
    }
  }
  if (par.empty()) {
    WarnAsErrorIfSet("No PARs available for chromosome {}", chrom);
  }
  if (name.empty()) {
    name = chrom;
    WarnAsErrorIfSet("Using default name for chromosome {}: {}", chrom, name);
  }
  return std::make_pair(name, par);
}

/**
 * @brief Extract sex chromsosome information, predicted sex, and chromosomal median DP values from a VCF file.
 * @param param CLI parameters for `filter_variants`
 * @return Tuple: predicted sex, chrX name, chrY name, chrX PAR intervals, chrY PAR intervals, chromosomal median DPs
 */
static std::tuple<Sex, std::string, std::string, vec<Interval>, vec<Interval>, StrUnorderedMap<u32>> GetSexInfo(
    const FilterVariantsParam& param) {
  // Extract the name and PARs for chromosomes X and Y
  auto [chr_x_name, chr_x_par] = ExtractSexChromPAR(param.par_bed_x, kDefaultChrXName);
  auto [chr_y_name, chr_y_par] = ExtractSexChromPAR(param.par_bed_y, kDefaultChrYName);
  bool determine_sex = !chr_x_par.empty() && !chr_x_name.empty();
  auto sex = Sex::kUnknown;

  Logging::Info("Extracting DP values from VCF file {}", param.vcf_file);
  auto chr_to_dp = GetChromosomeMedianDP(param.vcf_file);
  if (chr_to_dp.empty()) {
    error::Error("No median DP values extracted");
  } else {
    Logging::Info("Extracted median DP for {} chromosomes", chr_to_dp.size());
    vec<std::string> items;
    items.reserve(chr_to_dp.size());
    for (const auto& [chrom, dp] : chr_to_dp) {
      items.emplace_back(chrom + ":" + std::to_string(dp));
    }
    std::ranges::sort(items);

    // We can determine the biological sex of input sample based on the median DP ratio between chr1 and chrX.
    if (!chr_to_dp.contains(param.sd_chr_name) || chr_to_dp.at(param.sd_chr_name) == 0) {
      WarnAsErrorIfSet("Cannot extract median DP for {}", param.sd_chr_name);
      determine_sex = false;
    }
    if (!chr_to_dp.contains(chr_x_name) || chr_to_dp.at(chr_x_name) == 0) {
      WarnAsErrorIfSet("Cannot extract median DP for {}", chr_x_name);
      determine_sex = false;
    }
    if (determine_sex) {
      // A normal human diploid genome has 2 copies of chr1.
      // A biological female has 2 copies of chrX; its expected median DP ratio of chr1 to chrX is ~1.0.
      // A biological male has 1 copy of chrX and 1 copy of chrY; its expected median DP ratio of chr1 to chrX is ~2.0.
      const u32 chr_x_dp = chr_to_dp.at(chr_x_name);
      if (chr_x_dp > 0) {
        f64 chr_1_to_x_dp_ratio = chr_to_dp.at(param.sd_chr_name) / static_cast<f64>(chr_x_dp);
        Logging::Info("Median DP ratio between {} and {}: {}", param.sd_chr_name, chr_x_name, chr_1_to_x_dp_ratio);
        // default sex is unknown
        auto ratio = lround(chr_1_to_x_dp_ratio);
        if (ratio == 1) {
          Logging::Info("Sex: female");
          sex = Sex::kFemale;
        } else if (ratio == 2) {
          Logging::Info("Sex: male");
          sex = Sex::kMale;
        }
      }
    }
  }
  if (sex == Sex::kUnknown) {
    WarnAsErrorIfSet("Cannot determine the biological sex of input sample");
  }
  return std::make_tuple(sex, chr_x_name, chr_y_name, chr_x_par, chr_y_par, chr_to_dp);
}

#ifdef SOMATIC_ENABLE
// This function is designed to filter a given region of a GATK Mutect2 output VCF in a paired somatic tumor normal
// analysis. It is designed to use a single thread and be run in parallel for multiple regions simultaneously as part of
// a region based end to end filtering workflow. If there are variants that need to be filtered in the specified region
// of the genome, BAM and VCF features for the region are first computed before variants are read in and filtered using
// a LightGBM ML model to determine if they are indeed somatic variants. Variants that are deemed to not be somatic from
// this model are then passed to a second LightGBM ML model that attempts to determine if they represent germline
// variation between the sample and the reference genome or noise. Finally, all variants are written to the output
// buffer.
void FilterSomaticTumorNormalRegion(const GlobalContext& global_ctx,
                                    FlowContext& flow_ctx,
                                    WorkerContext& worker_ctx,
                                    const TargetRegion& region) {
  using enum VariantType;

  // Load values for worker and global contexts
  auto& reader = worker_ctx.filter_variants_reader;

  const auto& somatic_calculator = worker_ctx.calculators[0];
  const auto& germline_fail_calculator = worker_ctx.calculators[1];

  auto& out_records = flow_ctx.out_records;
  vec<io::VcfRecordPtr> fail_records;

  const auto& somatic_scoring_cols = global_ctx.model_config.scoring_cols;
  const auto& germline_scoring_cols = global_ctx.model_config.germline_fail_scoring_cols;
  const auto& hdr = global_ctx.hdr;
  const auto& bed_regions = global_ctx.bed_regions;
  const auto& interest_regions = global_ctx.interest_regions;
  const auto& normalize_targets = global_ctx.normalize_targets;
  const auto somatic_ml_snv_threshold = global_ctx.somatic_tn_snv_ml_threshold;
  const auto somatic_ml_indel_threshold = global_ctx.somatic_tn_indel_ml_threshold;
  const auto germline_fail_ml_threshold = global_ctx.germline_fail_ml_threshold;
  const auto tumor_support_threshold = global_ctx.tumor_support_threshold;

  if (!reader.SetRegion(region.chrom, region.start, region.end)) {
    return;
  }

  auto header_info = GetVcfHeaderInfo(hdr);
  if (!header_info.has_tumor_normal) {
    // Somatic TN Workflow but VCF doesn't have both tumor and normal
    return;
  }

  // Compute BAM and VCF Features here
  auto [vcf_features, bam_features, ref_features] =
      ComputeBamAndVcfFeaturesForRegion(global_ctx, worker_ctx, region, bed_regions, interest_regions);

  while (const auto& record = reader.GetNextRegionRecord(BCF_UN_ALL)) {
    const auto& chrom = record->Chromosome();
    if (record->Position() < 0) {
      WarnAsErrorIfSet("Found VCF record with position less than 0");
      continue;
    }
    const auto pos = static_cast<u64>(record->Position());

    if (chrom != region.chrom || pos < region.start || pos >= region.end) {
      continue;
    }

    std::optional<u32> normalize_target;
    if (!normalize_targets.empty()) {
      try {
        normalize_target = normalize_targets.at(chrom);
      } catch (const std::out_of_range& e) {
        // no normalize target for this chromosome
      }
    }

    const auto& ref = record->Allele(0);
    bool ref_acgt_only = !IsAnyNotACTG(ref);
    const auto& alt = record->Allele(1);
    auto record_copy = CopySomaticRecord(record, hdr);
    if (ref_acgt_only) {
      if (alt == "*") {
        out_records.emplace_back(record_copy);
      } else {
        const auto& [ref_trimmed, alt_trimmed] = TrimVariant(ref, alt);
        VariantId id(chrom, pos, ref_trimmed, alt_trimmed);
        VcfFeature vcf_feat;
        UnifiedVariantFeature var_feat;
        bool feat_found = false;
        try {
          vcf_feat = vcf_features.at(chrom).at(pos).at(id);
          var_feat = bam_features.at(chrom).at(pos).at(id);
          feat_found = vcf_feat.genotype != Genotype::kGTNA;
        } catch (std::out_of_range& e) {
          // VCF/BAM feature(s) not found for this variant
        }

        const UnifiedReferenceFeature& snv_ins_ref_feat = FindRefFeature(ref_features, chrom, pos);
        // del reference feature needs to be offset by 1
        const UnifiedReferenceFeature& del_ref_feat = FindRefFeature(ref_features, chrom, pos + 1);
        auto& ref_feat = id.type == kDeletion ? del_ref_feat : snv_ins_ref_feat;

        record_copy->SetAlleles({ref_trimmed, alt_trimmed});
        if (feat_found) {
          const auto& feature_vec =
              GetFeatureVec(somatic_scoring_cols, id, var_feat, ref_feat, vcf_feat, normalize_target);
          auto score = somatic_calculator.CalculateScore(feature_vec);
          SetSomaticTNValues(
              record_copy, id, var_feat, ref_feat, vcf_feat, header_info.tumor_index < header_info.normal_index, score);
          auto tumor_support = var_feat.tumor_support;
          auto ml_threshold_to_use = id.type == kSNV ? somatic_ml_snv_threshold : somatic_ml_indel_threshold;
          if (score < ml_threshold_to_use || tumor_support < tumor_support_threshold) {
            if (score < ml_threshold_to_use) {
              FailSomaticTNRecord(record_copy, kFilteringFailSomaticTNId);
            }
            if (tumor_support < tumor_support_threshold) {
              FailSomaticTNRecord(record_copy, kFilteringCountsId);
            }
            fail_records.emplace_back(record_copy);
          }
          out_records.emplace_back(record_copy);

        } else {
          // features not found; copy the record and set FILTER to FAIL
          FailSomaticTNRecord(record_copy, kFilteringFailId);
          out_records.emplace_back(record_copy);
        }
      }
    } else {
      // REF/ALT contains non-ACGT; copy the record and set FILTER to FAIL
      FailSomaticTNRecord(record_copy, kFilteringFailId);
      out_records.emplace_back(record_copy);
    }
  }

  // Once first filtering pass has been completed, iterate through the out_records again
  // This time apply germline-fail model only to the failing records
  for (const auto& record : fail_records) {
    const auto& chrom = record->Chromosome();
    const auto pos = static_cast<u64>(record->Position());

    const auto& ref = record->Allele(0);
    const auto& alt = record->Allele(1);
    // Know that these records have features based on how they were added to the vector above
    VariantId id(chrom, pos, ref, alt);
    VcfFeature vcf_feat;
    UnifiedVariantFeature var_feat;
    vcf_feat = vcf_features.at(chrom).at(pos).at(id);
    var_feat = bam_features.at(chrom).at(pos).at(id);

    const UnifiedReferenceFeature& snv_ins_ref_feat = FindRefFeature(ref_features, chrom, pos);
    // del reference feature needs to be offset by 1
    const UnifiedReferenceFeature& del_ref_feat = FindRefFeature(ref_features, chrom, pos + 1);
    auto& ref_feat = id.type == kDeletion ? del_ref_feat : snv_ins_ref_feat;

    std::optional<u32> normalize_target;
    if (!normalize_targets.empty()) {
      try {
        normalize_target = normalize_targets.at(chrom);
      } catch (const std::out_of_range& e) {
        // no normalize target for this chromosome
      }
    }

    const auto& feature_vec = GetFeatureVec(germline_scoring_cols, id, var_feat, ref_feat, vcf_feat, normalize_target);
    auto score = germline_fail_calculator.CalculateScore(feature_vec);
    if (score < germline_fail_ml_threshold) {
      FailSomaticTNRecord(record, kFilteringFailSomaticTNGermlineId);
    }
  }
}

// This function is designed to filter a given region of a GATK Mutect2 output VCF in a tumor only somatic analysis. It
// is designed to use a single thread and be run in parallel for multiple regions simultaneously as part of a region
// based end to end filtering workflow. If there are variants in the specified region of the genome, BAM and VCF
// features for the region are first computed before variants are read in and filtered using a LightGBM ML model to
// determine if they are indeed somatic variants. If variants are deemed to be phased by GATK Mutect2 a separate set of
// filtering criteria is applied.
void FilterSomaticRegion(const GlobalContext& global_ctx,
                         FlowContext& flow_ctx,
                         WorkerContext& worker_ctx,
                         const TargetRegion& region) {
  // Load values for worker and global contexts
  auto& reader = worker_ctx.filter_variants_reader;

  const auto& somatic_calculator = worker_ctx.calculators[0];

  auto& out_records = flow_ctx.out_records;

  const auto& scoring_cols = global_ctx.model_config.scoring_cols;
  const auto& hdr = global_ctx.hdr;
  const auto& bed_regions = global_ctx.bed_regions;
  const auto& interest_regions = global_ctx.interest_regions;
  const auto& normalize_targets = global_ctx.normalize_targets;
  const bool make_forcecalls = global_ctx.force_calls.has_value();
  const bool phased = global_ctx.phased;
  const auto& force_calls = global_ctx.force_calls.has_value() ? global_ctx.force_calls.value() : vec<BedRegion>{};
  const auto& block_list = global_ctx.block_list.has_value() ? global_ctx.block_list.value() : StrUnorderedSet{};
  const auto& hotspots = global_ctx.hotspots.has_value() ? global_ctx.hotspots.value() : StrUnorderedSet{};

  if (!reader.SetRegion(region.chrom, region.start, region.end)) {
    return;
  }
  // Compute BAM and VCF Features here
  auto [vcf_features, bam_features, ref_features] =
      ComputeBamAndVcfFeaturesForRegion(global_ctx, worker_ctx, region, bed_regions, interest_regions);

  // Compute thresholds and setup filtering settings for somatic filtering
  auto weight_counts_thresholds =
      CalculateWeightedCountThresholdsPerSubstitutionType(bam_features, 0, global_ctx.weighted_counts_threshold);
  const auto filter_settings = FilterSettings(global_ctx.bam_feat_params.min_mapq,
                                              global_ctx.bam_feat_params.min_bq,
                                              weight_counts_thresholds,
                                              global_ctx.ml_threshold,
                                              block_list,
                                              global_ctx.hotspot_weighted_counts_threshold,
                                              global_ctx.hotspot_ml_threshold,
                                              global_ctx.min_allele_freq_threshold,
                                              hotspots);
  const auto phased_filter_settings = PhasedFilterSettings(global_ctx.bam_feat_params.min_mapq,
                                                           global_ctx.bam_feat_params.min_bq,
                                                           block_list,
                                                           global_ctx.min_phased_allele_freq,
                                                           global_ctx.max_phased_allele_freq,
                                                           global_ctx.min_alt_counts);

  u64 prev_pos = region.start;

  StrUnorderedSet seen_forced_calls;
  while (const auto& record = reader.GetNextRegionRecord(BCF_UN_ALL)) {
    const auto chrom = record->Chromosome();
    if (record->Position() < 0) {
      WarnAsErrorIfSet("Found VCF record with position less than 0");
      continue;
    }
    const auto pos = static_cast<u64>(record->Position());
    const auto pos_int = static_cast<s32>(record->Position());
    if (chrom != region.chrom || pos < region.start || pos >= region.end) {
      continue;
    }

    auto ref = record->Allele(0);
    auto alt = record->Allele(1);
    if (record->NumAlleles() > 2) {
      throw error::Error(
          fmt::format("VCF file contains a variant with more than 2 alleles: {}:{}. "
                      "Run 'gatk LeftAlignAndTrimVariants -R <genome> -V <this VCF> -O <new VCF> "
                      "--split-multi-allelics' to normalize the VCF file.",
                      chrom,
                      pos + 1));
    }
    auto key = GetVariantCorrelationKey(chrom, pos, ref, alt);
    auto key_unpadded = GetVariantCorrelationKey(chrom, pos, ref, alt, false);
    // If the record has a matched VCF and BAM feature, get it here
    VariantId vid(chrom, pos, ref, alt);
    VcfFeature vcf_feature;
    if (!vcf_features.empty()) {
      try {
        vcf_feature = vcf_features.at(chrom).at(pos).at(vid);
      } catch (const std::out_of_range& e) {
      }
    }
    const u64 ref_feature_pos = VariantId::GetRefFeaturePos(ref, alt, pos);
    const UnifiedReferenceFeature& reference_feature = ref_features[chrom][ref_feature_pos];
    // Must have a matching BAM feature for somatic calling. VCF feature is not required but may be used if VCF features
    // have been specified
    UnifiedVariantFeature variant_feature;
    bool found_variant = false;
    if (bam_features[chrom][pos].contains(vid)) {
      variant_feature = bam_features[chrom][pos][vid];
      found_variant = true;
    }
    if (found_variant) {
      std::optional<u32> normalize_target;
      if (!normalize_targets.empty()) {
        try {
          normalize_target = normalize_targets.at(chrom);
        } catch (const std::out_of_range& e) {
          // no normalize target for this chromosome
        }
      }
      auto feature_vec =
          GetFeatureVec(scoring_cols, vid, variant_feature, reference_feature, vcf_feature, normalize_target);
      variant_feature.ml_score = somatic_calculator.CalculateScore(feature_vec);
      variant_feature.filter_status = FilterVariant(vid, variant_feature, reference_feature, filter_settings);
    }
    if (make_forcecalls) {
      // For forced calling we must check every position from the last seen variant position up to the current record's
      // position. Any force calls that need to be made should be done so before dealing with the current record to keep
      // the sorted order intact.
      u64 max_pos_to_check = pos;
      // TODO: Handle case where there are no VCF records for a chromosome but there are features and forced
      // sites Check if we have any remaining variants in force set
      for (u64 i = prev_pos; i < max_pos_to_check; i++) {
        auto i_int = static_cast<s32>(i);
        if (bam_features[chrom].contains(i) &&
            std::ranges::find(force_calls, BedRegion(chrom, i_int, i_int + 1)) != force_calls.end() &&
            !seen_forced_calls.contains(chrom + ":" + std::to_string(i))) {
          for (const auto& [_vid, _feat] : bam_features[chrom][i]) {
            io::VcfRecordPtr new_record = io::VcfRecord::CreateFromHeader(hdr);
            new_record->SetChromosome(chrom);
            new_record->SetPosition(i_int);
            SetupForceRecord(new_record, _vid);
            auto key_to_check = GetVariantCorrelationKey(chrom, i, _vid.ref, _vid.alt, false);
            auto is_hotspot = filter_settings.hotspots.contains(key_to_check);
            const UnifiedReferenceFeature& ref_feature = ref_features[chrom][i];
            UpdateRecord(new_record, _feat, ref_feature, is_hotspot);
            out_records.emplace_back(new_record);
          }
        }
      }
      prev_pos = pos;
    }
    bool phased_record = false;
    if (phased) {
      vec<std::string> pid_values;
      try {
        // If PID or PGT fields not in VCF record, a runtime error is thrown
        pid_values = record->GetFormatField<std::string>("PID");
      } catch (std::runtime_error& e) {
        // No field for phased variant calls in the header, Don't need to do anything
      }
      if (!pid_values.empty()) {
        phased_record = true;
      }
    }
    auto record_copy = record->Clone(hdr);
    bool found = false;
    if (!phased_record && found_variant) {
      for (const auto& item : variant_feature.filter_status) {
        record_copy->AddFilter(item);
      }
      auto is_hotspot = filter_settings.hotspots.contains(key_unpadded);
      UpdateRecord(record_copy, variant_feature, reference_feature, is_hotspot);
      out_records.emplace_back(record_copy);
      found = true;
    } else if (phased_record) {
      // Phased variant
      found = true;
      vec<s32> n_ad = record->GetFormatField<s32>("AD");
      vec<f32> n_af = record->GetFormatField<f32>("AF");
      vec<s32> n_mq = record->GetInfoField<s32>("MMQ");
      vec<s32> n_bq = record->GetInfoField<s32>("MBQ");
      vec<std::string> filters = FilterPhasedVariant(key,
                                                     phased_filter_settings,
                                                     static_cast<u32>(n_ad[1]),
                                                     n_af[0],
                                                     *std::ranges::max_element(n_mq),
                                                     *std::ranges::max_element(n_bq));
      for (const auto& item : filters) {
        record_copy->AddFilter(item);
      }
      UpdatePhasedRecord(record_copy);
      out_records.emplace_back(record_copy);
    }
    if (make_forcecalls && found &&
        std::ranges::find(force_calls, BedRegion(chrom, pos_int, pos_int + 1)) != force_calls.end()) {
      seen_forced_calls.insert(chrom + ":" + std::to_string(pos));
    }
  }

  // Check if there are any forcecalls left in the region after the last VCF variant.
  if (make_forcecalls) {
    auto max_pos_to_check = region.end;
    // TODO: Handle case where there are no VCF records for a chromosome but there are features and forced
    // sites Check if we have any remaining variants in force set
    for (u64 i = prev_pos; i < max_pos_to_check; i++) {
      auto i_int = static_cast<s32>(i);
      if (bam_features[region.chrom].contains(i) &&
          std::ranges::find(force_calls, BedRegion(region.chrom, i_int, i_int + 1)) != force_calls.end() &&
          seen_forced_calls.contains(region.chrom + ":" + std::to_string(i))) {
        for (const auto& [_vid, _feat] : bam_features[region.chrom][i]) {
          io::VcfRecordPtr new_record = io::VcfRecord::CreateFromHeader(hdr);
          new_record->SetChromosome(region.chrom);
          new_record->SetPosition(i_int);
          SetupForceRecord(new_record, _vid);
          auto key_to_check = GetVariantCorrelationKey(region.chrom, i, _vid.ref, _vid.alt, false);
          auto is_hotspot = filter_settings.hotspots.contains(key_to_check);
          const UnifiedReferenceFeature& ref_feature = ref_features[region.chrom][i];
          UpdateRecord(new_record, _feat, ref_feature, is_hotspot);
          out_records.emplace_back(new_record);
        }
      }
    }
  }
}
#endif  // SOMATIC_ENABLE

void FilterVariantsClass::SetGlobalContext() {
  // The GlobalContext struct contains various parameters required for workflow-specific filtering, and they do not
  // change between regions. The exact parameters may differ for each workflow, but certain parameters are only
  // applicable to specific workflows.

  const ComputeBamFeaturesParams bam_feat_params{
      .feature_cols = _model_config.feature_cols,
      .min_bq = _param.min_bq,
      .min_mapq = _param.min_mapq,
      .min_allowed_distance_from_end = _param.min_allowed_distance_from_end,
      .min_family_size = _min_family_size,
      .max_read_variant_count = _param.max_read_variant_count,
      .max_read_variant_count_normalized = _param.max_read_variant_count_normalized,
      .min_homopolymer_length = _param.min_homopolymer_length,
      .duplex = _param.duplex,
      .filter_homopolymer = _param.filter_homopolymer,
      .tumor_read_group = std::nullopt,
      .decode_yc = _param.decode_yc,
      .min_base_type = _param.min_base_type};

  switch (_model_config.workflow) {
    case Workflow::kGermlineMultiSample:
    case Workflow::kGermline: {
      auto [sex, chr_x_name, chr_y_name, chr_x_par, chr_y_par, chr_to_dp] = GetSexInfo(_param);
      _global_ctx.genome = _param.genome;
      _global_ctx.bam_feat_params = bam_feat_params;
      _global_ctx.ref_seqs = _ref_seqs;
      _global_ctx.bed_regions =
          _param.bed_file.has_value() ? GetChromIntervalMap(_param.bed_file).value() : ChromIntervalsMap{};
      _global_ctx.interest_regions = _param.interest_bed_file.has_value()
                                         ? GetChromIntervalMap(_param.interest_bed_file).value()
                                         : ChromIntervalsMap{};
      _global_ctx.model_config = _model_config;
      _global_ctx.hdr = _hdr;
      _global_ctx.normalize_targets = _param.normalize_features ? chr_to_dp : StrUnorderedMap<u32>();
      _global_ctx.sex = sex;
      _global_ctx.chr_x_name = chr_x_name;
      _global_ctx.chr_y_name = chr_y_name;
      _global_ctx.chr_x_par = chr_x_par;
      _global_ctx.chr_y_par = chr_y_par;
      _global_ctx.skip_variants_vcf = _param.skip_variants_vcf;
    }
#ifdef SOMATIC_ENABLE
    case Workflow::kSomatic: {
      StrUnorderedSet block_list;
      if (param.block_list) {
        block_list = LoadBlockList(*param.block_list);
      }
      StrUnorderedSet hotspots;
      if (param.hotspot_list) {
        hotspots = LoadHotspotVariants(*param.hotspot_list);
      }
      // Force Calls BED is 0-based. VCF positions and features in file are 1-based
      // When read in VCF Position and features are 0-based
      vec<BedRegion> force_calls;
      if (param.forcecall_list) {
        force_calls = io::ParseBedFile(*param.forcecall_list);
      }

      GlobalContext g_ctx{
          .genome = param.genome,
          .bam_feat_params = bam_feat_params,
          .ref_seqs = ref_seqs,
          .bed_regions = param.bed_file.has_value() ? GetChromIntervalMap(param.bed_file).value() : ChromIntervalsMap{},
          .interest_regions = param.interest_bed_file.has_value() ? GetChromIntervalMap(param.interest_bed_file).value()
                                                                  : ChromIntervalsMap{},
          .model_config = model_config,
          .hdr = hdr,
          .normalize_targets =
              param.normalize_features ? GetChromosomeMedianDP(param.vcf_file) : StrUnorderedMap<u32>(),
          .phased = param.phased,
          .force_calls = force_calls,
          .hotspots = hotspots,
          .block_list = block_list,
          .min_allele_freq_threshold = param.min_allele_freq_threshold,
          .weighted_counts_threshold = param.weighted_counts_threshold,
          .hotspot_weighted_counts_threshold = param.hotspot_weighted_counts_threshold,
          .ml_threshold = param.ml_threshold,
          .hotspot_ml_threshold = param.hotspot_ml_threshold,
          .min_phased_allele_freq = param.min_phased_allele_freq,
          .max_phased_allele_freq = param.max_phased_allele_freq,
          .min_alt_counts = param.min_alt_counts,
          .skip_variants_vcf = param.skip_variants_vcf};
      global_ctx = g_ctx;
    }
    case Workflow::kSomaticTumorNormal: {
      GlobalContext g_ctx{
          .genome = param.genome,
          .bam_feat_params = bam_feat_params,
          .ref_seqs = ref_seqs,
          .bed_regions = param.bed_file.has_value() ? GetChromIntervalMap(param.bed_file).value() : ChromIntervalsMap{},
          .model_config = model_config,
          .hdr = hdr,
          .normalize_targets =
              param.normalize_features ? GetChromosomeMedianDP(param.vcf_file) : StrUnorderedMap<u32>(),
          .somatic_tn_snv_ml_threshold = param.somatic_tn_snv_ml_threshold,
          .somatic_tn_indel_ml_threshold = param.somatic_tn_indel_ml_threshold,
          .tumor_support_threshold = param.tumor_support_threshold,
      };
      global_ctx = g_ctx;
    }
#endif  // SOMATIC_ENABLE
    default: {
      break;
    }
  }
}

void FilterVariantsClass::ParallelFiltering(const io::VcfWriter& out_file) {
  // Creates a global context and worker contexts for each worker thread and uses taskflow to execute parallel filtering
  // and output writing tasks, filtering regions of the input VCF in parallel and writing VCF records out in sorted
  // order.

  // Initialize the global context, these values are shared across all regions and all tasks, they are read-only.
  SetGlobalContext();
  if (auto num_warnings = CheckVcfFeatureResources(
          _model_config.vcf_feature_cols, _global_ctx.interest_regions, _hdr, _param.pop_af_vcf.has_value());
      num_warnings > 0) {
    WarnAsErrorIfSet(
        "There were {} warnings while checking VCF feature resources; VCF feature extraction may not be accurate!",
        num_warnings);
  }
  Progress progress(_partitioned_regions.size());
  // Allocate a window of flow contexts, a flow context is a temporary storage for output records for a region.
  // A window is used to help address uneven processing times for different regions, and to reduce the memory
  // footprint of the program. A larger window size means that more intermediate records are stored in memory,
  // but also that more processing can happen in parallel before stalling.
  const u32 flow_ctx_window_size = _param.filter_variant_window_size;
  vec<FlowContext> flow_ctxs{flow_ctx_window_size};

  tf::Taskflow taskflow;
  tf::Executor executor{_param.threads};

  vec<tf::Task> filter_tasks;
  filter_tasks.reserve(_partitioned_regions.size());

  vec<tf::Task> write_tasks;
  write_tasks.reserve(_partitioned_regions.size());

  // For each thread setup the required resources
  _workers.reserve(_param.threads);
  for (u32 i = 0; i < _param.threads; ++i) {
    try {
      auto worker_ctx = std::make_unique<WorkerContext>(
          _param.vcf_file, _param.pop_af_vcf, _param.bam_files, _model_config, _alignment_reader_cache);
      _workers.emplace_back(_global_ctx, worker_ctx, _hdr);
    } catch (std::exception& e) {
      throw error::Error("Error creating WorkerContext: {}", e.what());
    }
  }

  // For each region we will create a filter task and an output writing task, the filter task will filter
  // and create new VCF records stored in the flow context. The write task will write the output records to
  // the output VCF file.
  for (size_t i = 0; i < _partitioned_regions.size(); i++) {
    const auto& region = _partitioned_regions.at(i);
    // Determine the flow context to use for this region, this is done by cycling through the flow context
    // a flow context will be reused over and over.
    auto flow_ctx_idx = i % flow_ctx_window_size;
    auto& flow_ctx = flow_ctxs.at(flow_ctx_idx);

    // Create a filter task for the given region, a filter task is dependent on the write task from the previous
    // window which has the same flow_ctx_idx, this ensures that this filter task does not start until the previous
    // window has been written to the output file.
    auto filter_task = taskflow.emplace([this, &flow_ctx, &region, &executor, &progress]() {
      try {
        auto& worker = _workers.at(executor.this_worker_id());
        progress.UpdateAndLog(log::LogLevel::kInfo);
        // Filtering tasks will differ based on workflow
        switch (_model_config.workflow) {
          case Workflow::kGermlineMultiSample:
          case Workflow::kGermline: {
            worker.FilterGermlineRegion(region, flow_ctx);
            break;
          }
#ifdef SOMATIC_ENABLE
          case Workflow::kSomatic: {
            FilterSomaticRegion(global_ctx, flow_ctx, worker_ctx, region);
            break;
          }
          case Workflow::kSomaticTumorNormal: {
            FilterSomaticTumorNormalRegion(global_ctx, flow_ctx, worker_ctx, region);
            break;
          }
#endif  // SOMATIC_ENABLE
          default: {
            break;
          }
        }
      } catch (std::exception& e) {
        throw error::Error(
            "Error filtering variants in region '{}:{}-{}': {}", region.chrom, region.start, region.end, e.what());
      }
    });
    filter_tasks.emplace_back(filter_task);
    if (i >= flow_ctx_window_size) {
      filter_task.succeed(write_tasks.at(i - flow_ctx_window_size));
    }

    // Create a write task for the given region, a write task will have 2 dependencies:
    //   1. The filter task for this region, this ensures that the flow context has all the output records for this
    //   region.
    //   2. The previous write task, this ensures that the output records are written in order.
    auto write_task = taskflow.emplace([&flow_ctx, &out_file, &region]() {
      try {
        for (const auto& record : flow_ctx.out_records) {
          out_file.WriteRecord(record);
        }
        // Clear out the records to ensure that the next region does not have any records from this region.
        flow_ctx.out_records.clear();
      } catch (std::exception& e) {
        throw error::Error("Error writing region '{}:{}-{}': {}", region.chrom, region.start, region.end, e.what());
      }
    });
    if (!write_tasks.empty()) {
      write_task.succeed(write_tasks.back());
    }
    write_tasks.emplace_back(write_task);

    write_task.succeed(filter_task);
  }

  Logging::Info("Filtering and writing '{}' regions...", progress.total_region_count);
  executor.run(taskflow).get();

  if (_partitioned_regions.size() != progress.total_region_count) {
    throw error::Error("Number of regions {} and tasks executed {} are not the same!",
                       _partitioned_regions.size(),
                       std::to_string(progress.total_region_count));
  }
}

#ifdef SOMATIC_ENABLE
/**
 * @brief Entry point for somatic tumor normal filtering. Model features for the somatic tn and germline fail model are
 * validated, somatic tn specific header information is added to the output VCF header and input files are validated
 * before parallel filtering is performed
 * @param param Input parameters
 * @param model_config A SVCConfig struct
 */
void FilterSomaticTN(const FilterVariantsParam& param, const SVCConfig& model_config) {
  VerifyModelCompatibility(model_config.germline_fail_model_file, model_config.germline_fail_scoring_cols);
  VerifyModelCompatibility(model_config.model_file, model_config.scoring_cols);

  vec<TargetRegion> partitioned_regions = PartitionVcfRegions(param.vcf_file, param.threads);
  // load the VCF File
  Logging::Info("Loading VCF file {}", param.vcf_file);
  io::VcfReader reader(param.vcf_file);
  if (1 != reader.GetHeader()->GetNumSamples()) {
    throw error::Error("VCF file must contain exactly one sample");
  }

  Logging::Info("Writing output VCF file {}", param.vcf_output.string());

  // Add somatic-tn specific header information
  auto hdr = reader.GetHeader()->Clone();
  if (param.command_line) {
    hdr->AddCustomMetaDataLine(GetVcfHeaderLine(*param.command_line));
  }
  hdr->AddFilterLine({kFilteringPassId, kFilteringGermlinePassDesc});
  hdr->AddFilterLine({kFilteringFailId, kFilteringGermlineFailDesc});
  hdr->AddFilterLine({kFilteringFailSomaticTNId, kFilteringFailSomaticTNDesc});
  hdr->AddFilterLine({kFilteringFailSomaticTNGermlineId, kFilteringFailSomaticTNGermlineDesc});
  hdr->AddFilterLine({kFilteringCountsId, kFilteringCountsDesc});
  // TODO : Confirm which of the INFO and FORMAT fields are already present in the hdr for somatic TN
  hdr->AddInfoLine({kPredId, kPredDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddInfoLine({kRefBQId, kRefBQDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddInfoLine({kRefMQId, kRefMQDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddInfoLine({kAltBQId, kAltBQDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddInfoLine({kAltMQId, kAltMQDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddInfoLine({kSubtypeId, kSubtypeDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddInfoLine({kContextId, kContextDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddInfoLine({kFieldPopaf, kPopafDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddFormatLine({kBaseqQualId, kBaseQualDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddFormatLine({kMapQualId, kMapQualDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddFormatLine({kDistanceId, kDistanceDesc, io::kNumberOne, io::FieldType::kFloat});

  if (!param.vcf_output.parent_path().empty()) {
    fs::create_directories(param.vcf_output.parent_path());
  }
  io::VcfWriter out_file(param.vcf_output, hdr);
  out_file.WriteHeader();

  StrSet chrom_names{};
  for (const auto& r : partitioned_regions) {
    chrom_names.insert(r.chrom);
  }
  io::FastaReader fasta_reader(param.genome);
  StrMap<std::string> ref_seqs;
  for (const auto& [chrom, length] : hdr->GetContigLengths()) {
    if (chrom_names.contains(chrom)) {
      ref_seqs[chrom] = fasta_reader.GetSequence(chrom, 0, static_cast<s32>(length));
    }
  }

  if (param.bam_files.empty()) {
    throw error::Error("No BAM file provided");
  }

  u32 min_family_size = param.min_family_size;
  if (param.duplex && param.min_family_size > 2) {
    WarnAsErrorIfSet("Lowering min family size from {} to 2 for duplex data.", min_family_size);
    min_family_size = 2;
  }

  ParallelFiltering(param, model_config, partitioned_regions, hdr, min_family_size, ref_seqs, out_file);

  Logging::Info("Wrote output VCF file {}", param.vcf_output);
}
#endif  // SOMATIC_ENABLE

void FilterVariantsClass::SetReferenceSequences() {
  StrSet chrom_names{};
  for (const auto& r : _partitioned_regions) {
    // Only store the chromosomes that have variants in the VCF file that need to be processed
    chrom_names.insert(r.chrom);
  }
  io::FastaReader fasta_reader(_param.genome);
  for (const auto& [chrom, length] : _hdr->GetContigLengths()) {
    if (chrom_names.contains(chrom)) {
      _ref_seqs[chrom] = fasta_reader.GetSequence(chrom, 0, static_cast<s32>(length));
    }
  }
}

void FilterVariantsClass::FilterGermline() {
  // Verify that the model files and scoring columns are compatible. The specified scoring columns must match the
  // features used in the model both in name and order. Any mismatch means the models are not compatible with the
  // requested workflow based on the provided config and will not produce reliable results.
  VerifyModelCompatibility(_model_config.snv_model_file, _model_config.snv_scoring_cols);
  VerifyModelCompatibility(_model_config.indel_model_file, _model_config.indel_scoring_cols);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  if (1 != reader.GetHeader()->GetNumSamples()) {
    // Germline filtering is designed for a single sample only.
    throw error::Error("VCF file must contain exactly one sample");
  }

  if (!_param.vcf_output.parent_path().empty()) {
    auto parent_path = _param.vcf_output.parent_path();
    if (!fs::exists(parent_path)) {
      Logging::Info("Creating output directory {}", parent_path);
      fs::create_directories(parent_path);
    }
  }
  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());
  // Add germline specific header information
  _hdr = reader.GetHeader()->Clone();
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(GetVcfHeaderLine(*_param.command_line));
  }
  _hdr->AddFilterLine({kFilteringPassId, kFilteringGermlinePassDesc});
  _hdr->AddFilterLine({kFilteringFailId, kFilteringGermlineFailDesc});
  _hdr->AddFilterLine({kFilteringMissingFeatureId, kFilteringMissingFeatureDesc});
  _hdr->AddFilterLine({kFilteringFalsePositiveId, kFilteringFalsePositiveDesc});
  _hdr->AddFilterLine({kFilteringMultialleleFormatId, kFilteringMultialleleFormatDesc});
  _hdr->AddFilterLine({kFilteringMultiallelePartnerId, kFilteringMultiallelePartnerDesc});
  _hdr->AddFilterLine({kFilteringMultialleleConflictId, kFilteringMultialleleConflictDesc});
  _hdr->AddFilterLine({kFilteringNonAcgtRefAltId, kFilteringNonAcgtRefAltDesc});
  _hdr->AddFormatLine({kGermlineMLId, kGermlineMLDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kGermlineGATKDPId, kGermlineGATKDPDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kGermlineGATKGTId, kGermlineGATKGTDesc, io::kNumberOne, io::FieldType::kString});
  _hdr->AddFormatLine({kGermlineGATKADId, kGermlineGATKADDesc, io::kNumberDot, io::FieldType::kInteger});
  _hdr->AddFormatLine({kGermlineGATKAltId, kGermlineGATKAltDesc, io::kNumberDot, io::FieldType::kString});
  _hdr->AddFormatLine({kGermlineGnomadAFId, kGermlineGnomadAFDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineRefAvgMapqId, kGermlineRefAvgMapqDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineAltAvgMapqId, kGermlineAltAvgMapqDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineRefAvgDistId, kGermlineRefAvgDistDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineAltAvgDistId, kGermlineAltAvgDistDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineDensity100BPId, kGermlineDensity100BPDesc, io::kNumberOne, io::FieldType::kInteger});
  const io::VcfWriter out_file(_param.vcf_output, _hdr);
  out_file.WriteHeader();

  // Partition the VCF file into regions for parallel processing
  _partitioned_regions = PartitionVcfRegions(_param.vcf_file, _param.threads);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  FilterVariantsClass::SetReferenceSequences();

  // Filter the input VCF in parallel
  ParallelFiltering(out_file);

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

#ifdef SOMATIC_ENABLE
/**
 * @brief Entry point for somatic tumor only filtering. A single model file is validated, BED regions parsed, BAM files
 * validated and somatic specific header lines added to the output VCF before parallel filtering and writing is
 * performed
 * @param param Input parameters
 * @param model_config A SVCConfig struct
 * @param min_family_size the minimum family size to use
 */
void FilterSomatic(const FilterVariantsParam& param, const SVCConfig& model_config, const u32 min_family_size) {
  VerifyModelCompatibility(param.model[0], model_config.scoring_cols);

  // load the VCF File
  Logging::Info("Loading VCF file {}", param.vcf_file);
  io::VcfReader reader(param.vcf_file);
  if (1 != reader.GetHeader()->GetNumSamples()) {
    throw error::Error("VCF file must contain exactly one sample");
  }

  // Need to use the BAM Header and idx to break up regions for somatic to match compute-bam-features region splitting
  // Once region splitting is unified with germline the use of the BAM header and idx can be removed here
  HtsFilePtr bam_file(io::HtsOpenFormat(param.bam_files[0].c_str(), "r", nullptr));
  SamHdrPtr bam_hdr(io::SamHdrRead(bam_file.get()));
  HtsIdxPtr idx(io::SamIdxLoad(bam_file.get(), param.bam_files[0].c_str()));

  StrMap<vec<Interval>> bed_regions;
  if (param.bed_file.has_value()) {
    // Use a vector of Intervals per chromosome vs BedRegion structs for storing BED regions. Makes lookups easier for
    // getting which intervals overlap a given position on a chromosomes
    // If there are bed regions, convert the vector of structs into a map so we can easily query
    auto regions = GetBedRegions(param.bed_file);
    if (regions.has_value()) {
      for (auto& region : regions.value()) {
        bed_regions[region.chromosome].emplace_back(region.start, region.end);
      }
      for (auto& [chrom, intervals] : bed_regions) {
        std::sort(intervals.begin(), intervals.end());
        bed_regions[chrom] = MergeIntervals(intervals);
      }
    }
  }
  // Regions that will be used to parallelize over are computed based on the BAM file to match compute_bam_features
  // This can be changed in the future
  vec<TargetRegion> partitioned_regions =
      PartitionRegionsForSomatic(param.max_bam_region_size_per_thread, bam_file, bam_hdr, idx, bed_regions);

  auto hdr = reader.GetHeader()->Clone();
  if (param.command_line) {
    hdr->AddCustomMetaDataLine(GetVcfHeaderLine(*param.command_line));
  }
  // Append descriptions for the FILTER fields
  hdr->AddFilterLine({kFilteringPassId, kFilteringPassDesc});
  hdr->AddFilterLine({kFilteringMapQualityId, kFilteringMapQualityDesc});
  hdr->AddFilterLine({kFilteringBlocklistedId, kFilteringBlocklistedDesc});
  hdr->AddFilterLine({kFilteringMLScoreId, kFilteringMLScoreDesc});
  hdr->AddFilterLine({kFilteringCountsId, kFilteringCountsDesc});
  hdr->AddFilterLine({kFilteringBaseQualityId, kFilteringBaseQualityDesc});
  hdr->AddFilterLine({kFilteringForcedId, kFilteringForcedDesc});
  hdr->AddFilterLine({kFilteringAFId, kFilteringAFDesc});
  hdr->AddFilterLine({kFilteringMinAltCountsId, kFilteringMinAltCountsDesc});

  // Append descriptions for the FORMAT fields
  hdr->AddFormatLine({kWeightedCountsId, kWeightedCountsDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  hdr->AddFormatLine({kMachineLearningId, kMachineLearningDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  hdr->AddFormatLine({kNonDuplexCountsId, kNonDuplexCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  hdr->AddFormatLine({kDuplexCountsId, kDuplexCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  hdr->AddFormatLine({kStrandBiasMetricId, kStrandBiasMetricDesc, io::kNumberOne, io::FieldType::kFloat});
  hdr->AddFormatLine({kPlusOnlyCountsId, kPlusOnlyCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  hdr->AddFormatLine({kMinusOnlyCountsId, kMinusOnlyCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  hdr->AddFormatLine({kSequenceContextId, kSequenceContextDesc, io::kNumberOne, io::FieldType::kString});
  hdr->AddFormatLine({kPhysicalPhasedId, kPhysicalPhasedDesc, io::kNumberOne, io::FieldType::kString});
  hdr->AddFormatLine({kHotspotId, kHotspotDesc, io::kNumberOne, io::FieldType::kInteger});

  StrSet chrom_names{};
  for (const auto& r : partitioned_regions) {
    chrom_names.insert(r.chrom);
  }
  io::FastaReader fasta_reader(param.genome);
  StrMap<std::string> ref_seqs;
  for (const auto& [chrom, length] : hdr->GetContigLengths()) {
    if (chrom_names.contains(chrom)) {
      ref_seqs[chrom] = fasta_reader.GetSequence(chrom, 0, static_cast<s32>(length));
    }
  }

  Logging::Info("Writing output VCF file {}", param.vcf_output);
  if (!param.vcf_output.parent_path().empty()) {
    fs::create_directories(param.vcf_output.parent_path());
  }
  io::VcfWriter out_file(param.vcf_output, hdr);
  out_file.WriteHeader();

  // Filter E2E by region using taskflow
  ParallelFiltering(param, model_config, partitioned_regions, hdr, min_family_size, ref_seqs, out_file);
}
#endif  // SOMATIC_ENABLE

void FilterVariantsClass::VerifyParameters() {
  // Model files must be specified from the CLI because `--model` is a required CLI parameter in `filter_variants`.
  // Therefore, do NOT retrieve model file paths from the config here.
  // The config contains the default output model file names intended for `train_model` only.
  if (_model_config.workflow == Workflow::kGermline || _model_config.workflow == Workflow::kGermlineMultiSample) {
    if (_param.model.size() == kGermlineNumModels) {
      _model_config.GetGermlineModelPaths(_param.model);
    } else {
      throw error::Error("Please specify two model files (SNV and indel) for the germline workflow");
    }
#ifdef SOMATIC_ENABLE
  } else if (model_config.workflow == Workflow::kSomatic) {
    // If using the somatic workflow there should only be one model passed in
    if (param.model.size() == 1) {
      model_config.model_file = param.model[0];
    } else {
      throw error::Error("Please specify one model file for the somatic workflow");
    }
  } else if (model_config.workflow == Workflow::kSomaticTumorNormal) {
    if (param.model.size() == 2) {
      model_config.GetSomaticTNModelPaths(param.model);
    } else {
      throw error::Error(
          "Please specify two model files (germline-fail and somatic) for the somatic-tumor-normal workflow");
    }
#endif  // SOMATIC_ENABLE
  }

  if (_param.bam_files.empty()) {
    throw error::Error("No BAM file provided");
  }
}

/**
 * The main entry point for the filter_variants module. This function performs checks to ensure the correct number of
 * models are passed based on the specified workflow before calling a workflow specific filtering function.
 * @param param A FilterVariantsParam struct that contains required input/output filepaths and filtering settings
 */
void FilterVariants(const FilterVariantsParam& param) {
  // This function acts as the entry point for the filter-variants tool. Based on the workflow the respective filtering
  // process is run

  // Create Filtering object
  FilterVariantsClass filtering(param);

  // Verify that the correct number of models have been passed based on the workflow, and that BAM input has been passed
  filtering.VerifyParameters();

  switch (param.workflow) {
    case Workflow::kGermlineMultiSample:
    case Workflow::kGermline: {
      filtering.FilterGermline();
      break;
    }
#ifdef SOMATIC_ENABLE
    case Workflow::kSomatic: {
      FilterSomatic(param, model_config, min_family_size);
      break;
    }
    case Workflow::kSomaticTumorNormal: {
      FilterSomaticTN(param, model_config);
      break;
    }
#endif  // SOMATIC_ENABLE
    default:
      throw error::Error("Specified Workflow has not yet been implemented");
  }
  Logging::Info("Filtered all variants");
}

}  // namespace xoos::svc
