#pragma once

#include <optional>
#include <string>
#include <vector>

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "core/command-line-info.h"
#include "core/config.h"
#include "core/genotype.h"
#include "core/workflow.h"
#include "core/yc-decode-method.h"
#include "filter-region.h"
#include "util/region-util.h"

namespace xoos::svc {

// CLI parameters for `filter_variants`
struct FilterVariantsParam {
  std::vector<fs::path> bam_files{};
  fs::path output_dir{};
  fs::path vcf_output{};
  fs::path genome{};
  fs::path model{};
  fs::path snv_model{};
  fs::path indel_model{};
  SVCConfig config{};
  std::optional<fs::path> config_file{};
  fs::path vcf_file{};
  std::optional<fs::path> bed_file{};
  std::optional<fs::path> block_list{};
  std::optional<fs::path> hotspot_list{};
  std::optional<fs::path> forcecall_list{};
  std::optional<fs::path> pop_af_vcf{};
  std::optional<ChromIntervalsMap> interest_regions{};
  u16 min_mapq{};
  u32 min_alt_counts{};
  f32 min_allele_freq_threshold{};
  f32 min_phased_allele_freq{};
  f32 max_phased_allele_freq{};
  u16 min_bq{};
  f32 min_allowed_distance_from_end{};
  u32 min_family_size{};
  u32 max_read_variant_count{};
  f32 max_read_variant_count_normalized{};
  size_t threads{};
  Workflow workflow{};
  f32 weighted_counts_threshold{};
  f32 hotspot_weighted_counts_threshold{};
  f32 ml_threshold{};
  std::optional<std::string> tumor_sample_name{};
  f32 somatic_tn_snv_ml_threshold{};
  f32 somatic_tn_indel_ml_threshold{};
  u32 tumor_support_threshold{};
  f32 hotspot_ml_threshold{};
  bool phased{};
  SequencingProtocol sequencing_protocol{SequencingProtocol::kDuplex};
  HomopolymerFilter filter_homopolymer{HomopolymerFilter::kNone};
  u32 min_homopolymer_length{};
  std::optional<CommandLineInfo> command_line;
  u32 max_bam_region_size_per_thread{};
  u32 max_vcf_region_size_per_thread{};
  u32 output_vcf_buffer_size{};
  FeatureNormalization normalize_features{FeatureNormalization::kNone};
  std::string sd_chr_name{};
  std::optional<fs::path> par_bed_x{};
  std::optional<fs::path> par_bed_y{};
  std::optional<fs::path> skip_variants_vcf{};
  YcDecodeMethod decode_yc{YcDecodeMethod::kNone};
  yc_decode::BaseType min_base_type;
};

const std::string kDefaultChrXName{"chrX"};
const std::string kDefaultChrYName{"chrY"};

const u32 kGermlineNumModels{2};
const u32 kDuplexMinFamilySize{2};

using GermlineScoreCalculator = std::function<Genotype(const vec<f64>&)>;

void FilterVariants(const FilterVariantsParam& param);

class FilterVariantsClass {
 private:
  FilterVariantsParam _param{};
  SVCConfig _model_config;
  u32 _min_family_size;
  io::VcfHeaderPtr _hdr;
  StrMap<std::string> _ref_seqs;
  vec<TargetRegion> _partitioned_regions;
  GlobalContext _global_ctx;
  vec<FilterRegionClass> _workers;
  AlignmentReaderCache _alignment_reader_cache;

 public:
  explicit FilterVariantsClass(const FilterVariantsParam& param)
      : _param(param), _model_config(_param.config), _min_family_size(param.min_family_size) {
    if (IsDuplexProtocol(param.sequencing_protocol) && param.min_family_size > kDuplexMinFamilySize) {
      _min_family_size = kDuplexMinFamilySize;
    }
  }

  /**
   * @brief Entry point for Germline filtering.
   *
   * Model files for SNV and Indel filtering are validated, before germline
   * specific VCF header lines are written to the output file, input BAMs are validated and parallel filtering is
   * performed.
   *
   * @note This function assumes that parameters have been verified by `VerifyParameters()` before being called.
   */
  void FilterGermline();

  void FilterGermlineTagging();

  /**
   * @brief Entry point for tumor only te filtering.
   *
   * A single model file is validated, BED regions parsed, BAM files
   * validated and somatic specific header lines added to the output VCF before parallel filtering and output VCF
   * writing is performed
   *
   * @note This function assumes that parameters have been verified by `VerifyParameters()` before being called.
   */
  void FilterTumorOnlyTe();

  /**
   * @brief Entry point for tumor normal filtering.
   *
   * A single model file is validated, BED regions parsed, BAM files
   * validated and somatic specific header lines added to the output VCF before parallel filtering and output VCF
   * writing is performed. NOTE: Currently only a single somatic model is used. Future implementations may also call
   * germline-tagging if the germline-tagging model is specified.
   *
   * @note This function assumes that parameters have been verified by `VerifyParameters()` before being called.
   */
  void FilterTumorNormal();

  /**
   * @brief Confirms that the tumor sample name is specified.
   */
  void CheckTumorNormalSampleName() const;

  /**
   * @brief Verifies input parameters required for the specified workflow.
   * @details Input parameters verified include:
   * - workflow is supported
   * - model files specified and compatible with the workflow
   * - at least one input BAM file
   * - tumor sample name specified for "germline-tagging", "tumor-only-te", and "tumor-normal-wgs" workflows
   */
  void VerifyParameters() const;

  /**
   * @brief Set up reference sequences for the chromosomes present in the partitioned regions. Only chromosomes where
   * there are variants to be filtered have their sequences pulled from the reference genome and stored in memory for
   * quick access.
   */
  void SetReferenceSequences();

  /**
   * @brief Set up the global context required for filtering variants. This sets parameters required for BAM feature
   * extraction and workflow specific filtering parameters that do not change between regions.
   *
   * @note This function assumes that _partitioned_regions has been set prior to being caled
   * @note This function assumes that _genome and _hdr have been set prior to being called and are well formed and not
   * null.
   */
  void SetGlobalContext();

  /**
   * @brief Perform parallel filtering of variants in the input VCF file based on the specified parameters and model
   * configuration. This function creates a global context and worker contexts for each thread, and uses taskflow to
   * execute parallel filtering and output writing tasks, filtering regions of the input VCF in parallel and writing VCF
   * records out in sorted order to the specified output file.
   * @param out_file A VCF output file to write records to.
   *
   * @note This function assumes that _partitioned_regions has been set and _ref_seqs populated prior to being called.
   */
  void ParallelFiltering(const io::VcfWriter& out_file);
};

}  // namespace xoos::svc
