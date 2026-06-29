#include "filter-variants.h"

#include <map>
#include <memory>
#include <utility>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/sex_predict/vcf-sex-predictor.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/container-functions.h>
#include <xoos/util/string-functions.h>

#include "compute-bam-features/compute-bam-features.h"
#include "compute-bam-features/progress-meter.h"
#include "compute-vcf-features/compute-vcf-features.h"
#include "compute-vcf-features/vcf-header-util.h"
#include "core/command-line-info.h"
#include "core/filtering.h"
#include "core/model-metadata.h"
#include "core/vcf-fields.h"
#include "util/file-util.h"
#include "util/log-util.h"
#include "util/parallel-compute-utils.h"

namespace xoos::svc {

using BedRegion = io::BedRegion;

/**
 * Loads a block list from file
 * @param block_list The blocklist in chr_pos_ref_alt format
 * @return The parsed left-padded blocklist
 */
static StrUnorderedSet LoadBlockList(const fs::path& block_list) {
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
 * @brief Helper function to extract sex chromosome name and PARs from a BED file.
 * @param bed_path Path of BED file containing PARs
 * @param default_chrom Default chromosome name if not specified in BED file
 * @return Chromosome name and vector of PARs.
 */
static std::pair<std::string, vec<Interval>> ExtractSexChromPAR(const std::optional<fs::path>& bed_path,
                                                                const std::string& default_chrom) {
  std::string name;
  vec<Interval> par;
  if (bed_path.has_value()) {
    auto chr_to_intervals = GetChromIntervalMap(bed_path).value();
    if (chr_to_intervals.empty()) {
      WarnAsErrorIfSet("No intervals found in {}", bed_path.value());
    } else if (chr_to_intervals.size() == 1) {
      name = chr_to_intervals.begin()->first;
      par = chr_to_intervals.begin()->second;
    } else {
      vec<std::string> ref_names;
      ref_names.reserve(chr_to_intervals.size());
      for (const auto& [chrom, intervals] : chr_to_intervals) {
        ref_names.emplace_back(chrom);
      }
      WarnAsErrorIfSet("PARs BED file has multiple reference names: {}", string::Join(ref_names, ","));
    }
  }
  if (par.empty()) {
    WarnAsErrorIfSet("No PARs available for chromosome {}", default_chrom);
  }
  if (name.empty()) {
    name = default_chrom;
    WarnAsErrorIfSet("Using default name for chromosome {}: {}", default_chrom, name);
  }
  return std::make_pair(name, par);
}

/**
 * @brief Extract sex chromsosome information, predicted sex, and chromosomal median DP values from a VCF file.
 * @param param CLI parameters for `filter_variants`
 * @return Tuple: predicted sex, chrX name, chrY name, chrX PAR intervals, chrY PAR intervals, chromosomal median DPs
 */
static std::tuple<sex_predict::Sex, std::string, std::string, vec<Interval>, vec<Interval>, ChromMedianDepth>
GetSexInfo(const FilterVariantsParam& param) {
  using enum sex_predict::Sex;
  // Extract the name and PARs for chromosomes X and Y
  auto [chr_x_name, chr_x_par] = ExtractSexChromPAR(param.par_bed_x, kDefaultChrXName);
  auto [chr_y_name, chr_y_par] = ExtractSexChromPAR(param.par_bed_y, kDefaultChrYName);
  auto sex = kUnknown;

  Logging::Info("Extracting DP values from VCF file {}", param.vcf_file);
  auto chr_med_dp = GetChromosomeMedianDepth(param.vcf_file);
  // use the "normal" sample's DP values for sex prediction
  if (chr_med_dp.normal.empty()) {
    error::Error("No median DP values extracted");
  }
  if (!chr_x_par.empty() && !chr_x_name.empty()) {
    sex = sex_predict::PredictSex(chr_med_dp.normal, param.sd_chr_name, chr_x_name);
    Logging::Info("Sex: {}", GetDescriptionForSex(sex));
    if (sex == kUnknown) {
      WarnAsErrorIfSet("Cannot determine the biological sex of input sample");
    }
  }
  return std::make_tuple(sex, chr_x_name, chr_y_name, chr_x_par, chr_y_par, chr_med_dp);
}

void FilterVariantsClass::SetGlobalContext() {
  using enum Workflow;
  // The GlobalContext struct contains various parameters required for workflow-specific filtering, and they do not
  // change between regions. The exact parameters may differ for each workflow, but certain parameters are only
  // applicable to specific workflows.

  _global_ctx.model_config = _model_config;

  // BAM feature extraction parameters
  _global_ctx.bam_feat_params =
      ComputeBamFeaturesParams{.feature_cols = _model_config.feature_cols,
                               .min_bq = _param.min_bq,
                               .min_mapq = _param.min_mapq,
                               .min_allowed_distance_from_end = _param.min_allowed_distance_from_end,
                               .min_family_size = _min_family_size,
                               .max_read_variant_count = _param.max_read_variant_count,
                               .max_read_variant_count_normalized = _param.max_read_variant_count_normalized,
                               .min_homopolymer_length = _param.min_homopolymer_length,
                               .sequencing_protocol = _param.sequencing_protocol,
                               .filter_homopolymer = _param.filter_homopolymer,
                               .tumor_sample_name = _param.tumor_sample_name,
                               .tumor_rg_ids = GetReadGroupIdsForSample(
                                   _alignment_reader_cache.Open(_param.bam_files, "r"), _param.tumor_sample_name),
                               .decode_yc = _param.decode_yc,
                               .min_base_type = _param.min_base_type};
  _global_ctx.skip_variants_vcf = _param.skip_variants_vcf;

  // reference genome and region parameters
  _global_ctx.genome = _param.genome;
  _global_ctx.ref_seqs = _ref_seqs;
  _global_ctx.bed_regions =
      _param.bed_file.has_value() ? GetChromIntervalMap(_param.bed_file).value() : ChromIntervalsMap{};
  _global_ctx.interest_regions =
      _param.interest_regions.has_value() ? _param.interest_regions.value() : ChromIntervalsMap{};

  _global_ctx.hdr = _hdr;
  std::tie(_global_ctx.vcf_info_metadata, _global_ctx.vcf_fmt_metadata) = _hdr->GetFieldMetadata();

  // Sex chromosome and normalization parameters
  const auto& [sex, chr_x_name, chr_y_name, chr_x_par, chr_y_par, chr_med_dp] = GetSexInfo(_param);
  if (_param.normalize_features == FeatureNormalization::kMedianDp) {
    // Check whether the calculated median depth values are valid for normalization.
    ValidateChromMedianDepthForNormalization(chr_med_dp, FeatureNamesHaveSampleContext(_model_config.workflow));
  }
  _global_ctx.normalize_targets =
      _param.normalize_features == FeatureNormalization::kMedianDp ? chr_med_dp : ChromMedianDepth{};
  _global_ctx.sex = sex;
  _global_ctx.chr_x_name = chr_x_name;
  _global_ctx.chr_y_name = chr_y_name;
  _global_ctx.chr_x_par = chr_x_par;
  _global_ctx.chr_y_par = chr_y_par;

  // Workflow specific parameters
  switch (_model_config.workflow) {
    case kGermlineMultiSample:
    case kGermline: {
      break;
    }
    case kTumorOnlyTe: {
      StrUnorderedSet block_list;
      if (_param.block_list) {
        block_list = LoadBlockList(*_param.block_list);
      }
      StrUnorderedSet hotspots;
      if (_param.hotspot_list) {
        hotspots = LoadHotspotVariants(*_param.hotspot_list);
      }
      _global_ctx.phased = _param.phased;
      _global_ctx.force_calls = GetChromIntervalMap(_param.forcecall_list);
      _global_ctx.hotspots = hotspots;
      _global_ctx.block_list = block_list;
      _global_ctx.min_allele_freq_threshold = _param.min_allele_freq_threshold;
      _global_ctx.weighted_counts_threshold = _param.weighted_counts_threshold;
      _global_ctx.hotspot_weighted_counts_threshold = _param.hotspot_weighted_counts_threshold;
      _global_ctx.ml_threshold = _param.ml_threshold;
      _global_ctx.hotspot_ml_threshold = _param.hotspot_ml_threshold;
      _global_ctx.min_phased_allele_freq = _param.min_phased_allele_freq;
      _global_ctx.max_phased_allele_freq = _param.max_phased_allele_freq;
      _global_ctx.min_alt_counts = _param.min_alt_counts;
      break;
    }
    case kGermlineTagging: {
      const auto tn_sample_idx = GetTumorNormalSampleIndexes(_hdr);
      if (tn_sample_idx.has_value()) {
        _global_ctx.vcf_tumor_index = tn_sample_idx->tumor_sample_idx;
        _global_ctx.vcf_normal_index = tn_sample_idx->normal_sample_idx;
      }
      _global_ctx.somatic_tn_snv_ml_threshold = _param.somatic_tn_snv_ml_threshold;
      _global_ctx.somatic_tn_indel_ml_threshold = _param.somatic_tn_indel_ml_threshold;
      _global_ctx.is_germline_tagging = true;
      break;
    }
    case kTumorNormalWgs: {
      _global_ctx.somatic_tn_snv_ml_threshold = _param.somatic_tn_snv_ml_threshold;
      _global_ctx.somatic_tn_indel_ml_threshold = _param.somatic_tn_indel_ml_threshold;
      _global_ctx.tumor_support_threshold = _param.tumor_support_threshold;
      const auto tn_sample_idx = GetTumorNormalSampleIndexes(_hdr);
      if (tn_sample_idx.has_value()) {
        _global_ctx.vcf_tumor_index = tn_sample_idx->tumor_sample_idx;
        _global_ctx.vcf_normal_index = tn_sample_idx->normal_sample_idx;
      }
      break;
    }
    default: {
      break;
    }
  }
}

/**
 * @brief Helper function to create score calculators based on workflow and model configuration.
 * @param param CLI parameters for `filter_variants`
 * @param model_config Model configuration for scoring
 * @return Vector of ScoreCalculator instances
 */
static vec<ScoreCalculator> CreateScoreCalculators(const FilterVariantsParam& param, const SVCConfig& model_config) {
  using enum Workflow;
  vec<ScoreCalculator> calculators;
  if (param.workflow == kGermline || param.workflow == kGermlineMultiSample) {
    calculators.emplace_back(
        param.snv_model, model_config.snv_scoring_cols.size(), model_config.snv_model_lgbm_prediction_params);
    calculators.emplace_back(
        param.indel_model, model_config.indel_scoring_cols.size(), model_config.indel_model_lgbm_prediction_params);
  } else {
    calculators.emplace_back(param.model, model_config.scoring_cols.size(), model_config.model_lgbm_prediction_params);
  }
  return calculators;
}

void FilterVariantsClass::ParallelFiltering(const io::VcfWriter& out_file) {
  using enum Workflow;
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
  const u32 flow_ctx_window_size = _param.output_vcf_buffer_size;
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
          _param.vcf_file, _param.pop_af_vcf, _param.bam_files, _alignment_reader_cache);
      worker_ctx->calculators = CreateScoreCalculators(_param, _model_config);
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
          case kGermlineMultiSample:
          case kGermline: {
            worker.FilterGermlineRegion(region, flow_ctx);
            break;
          }
          case kGermlineTagging: {
            worker.FilterGermlineTaggingRegion(region, flow_ctx);
            break;
          }
          case kTumorOnlyTe: {
            worker.FilterTumorOnlyTeRegion(region, flow_ctx);
            break;
          }
          case kTumorNormalWgs: {
            worker.FilterTumorNormalRegion(region, flow_ctx);
            break;
          }
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
  VerifyModelCompatibility(_param.snv_model, _model_config.snv_scoring_cols);
  VerifyModelCompatibility(_param.indel_model, _model_config.indel_scoring_cols);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  if (1 != reader.GetHeader()->GetNumSamples()) {
    // Germline filtering is designed for a single sample only.
    throw error::Error("VCF file must contain exactly one sample");
  }

  CreateParentDirectoryIfNotExists(_param.vcf_output);

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
  _hdr->AddFormatLine({kMachineLearningId, kMachineLearningDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineGATKDPId, kGermlineGATKDPDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kGermlineGATKGTId, kGermlineGATKGTDesc, io::kNumberOne, io::FieldType::kString});
  _hdr->AddFormatLine({kGermlineGATKADId, kGermlineGATKADDesc, io::kNumberDot, io::FieldType::kInteger});
  _hdr->AddFormatLine({kGermlineGATKAltId, kGermlineGATKAltDesc, io::kNumberDot, io::FieldType::kString});
  _hdr->AddFormatLine({kGermlineGnomadAFId, kGermlineGnomadAFDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineRefAvgMapqId, kGermlineRefAvgMapqDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineAltAvgMapqId, kGermlineAltAvgMapqDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineRefAvgDistId, kGermlineRefAvgDistDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineAltAvgDistId, kGermlineAltAvgDistDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kGermlineDensity100BPId, kGermlineDensity100BPDesc, io::kNumberOne, io::FieldType::kInteger});
  const io::VcfWriter out_file(_param.vcf_output, _hdr);
  out_file.WriteHeader();

  // Partition the VCF file into regions for parallel processing
  _partitioned_regions = PartitionVcfRegions(_param.vcf_file, _param.threads, std::nullopt);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  // Filter the input VCF in parallel
  ParallelFiltering(out_file);

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::FilterGermlineTagging() {
  VerifyModelCompatibility(_param.model, _model_config.scoring_cols);

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());
  // Add germline specific header information
  _hdr = reader.GetHeader()->Clone();
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(GetVcfHeaderLine(*_param.command_line));
  }
  _hdr->AddFormatLine({kMachineLearningId, kMachineLearningDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine(
      {kGermlineTaggingInfoGermlineId, kGermlineTaggingInfoGermlineDesc, io::kNumberZero, io::FieldType::kFlag});
  _hdr->AddInfoLine(
      {kGermlineTaggingInfoSomaticId, kGermlineTaggingInfoSomaticDesc, io::kNumberZero, io::FieldType::kFlag});
  if (!_hdr->HasInfoField(kFilteringPassId)) {
    _hdr->AddFilterLine({kFilteringPassId, kFilteringPassDesc});
  }
  if (!_hdr->HasInfoField(kFilteringFailId)) {
    _hdr->AddFilterLine({kFilteringFailId, kFilteringFailDesc});
  }
  const io::VcfWriter out_file(_param.vcf_output, _hdr);
  out_file.WriteHeader();

  // Partition the VCF file into regions for parallel processing
  _partitioned_regions = PartitionVcfRegions(_param.vcf_file, _param.threads, std::nullopt);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  // Filter the input VCF in parallel
  ParallelFiltering(out_file);

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::FilterTumorNormal() {
  VerifyModelCompatibility(_param.model, _model_config.scoring_cols);

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Loading VCF file {}", _param.vcf_file);
  const io::VcfReader reader(_param.vcf_file);
  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());
  // Add germline specific header information
  _hdr = reader.GetHeader()->Clone();
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(GetVcfHeaderLine(*_param.command_line));
  }
  // TODO: Replace these lines with a function that adds the appropriate headers for the tumor-normal workflow
  _hdr->AddFilterLine({kFilteringPassId, kFilteringGermlinePassDesc});
  _hdr->AddFilterLine({kFilteringFailId, kFilteringGermlineFailDesc});
  _hdr->AddFilterLine({kFilteringAFId, kFilteringAFDesc});
  _hdr->AddFilterLine({kFilteringFailSomaticTNMLId, kFilteringFailSomaticTNMLDesc});
  _hdr->AddFilterLine({kFilteringFailSomaticTNGermlineId, kFilteringFailSomaticTNGermlineDesc});
  _hdr->AddFilterLine({kFilteringCountsId, kFilteringCountsDesc});
  _hdr->AddFilterLine({kFilteringFailSomaticTNNormalADId, kFilteringFailSomaticTNNormalADDesc});
  _hdr->AddFilterLine({kFilteringMissingFeatureId, kFilteringMissingFeatureDesc});
  // TODO : Confirm which of the INFO and FORMAT fields are already present in the hdr for somatic TN
  _hdr->AddInfoLine({kRefBQId, kRefBQDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kRefMQId, kRefMQDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kAltBQId, kAltBQDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kAltMQId, kAltMQDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kSubtypeId, kSubtypeDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kContextId, kContextDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddInfoLine({kFieldPopaf, kPopafDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kMachineLearningId, kMachineLearningDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kBaseqQualId, kBaseQualDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kMapQualId, kMapQualDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kDistanceId, kDistanceDesc, io::kNumberOne, io::FieldType::kFloat});
  const io::VcfWriter out_file(_param.vcf_output, _hdr);
  out_file.WriteHeader();

  // Partition the VCF file into regions for parallel processing
  _partitioned_regions = PartitionVcfRegions(_param.vcf_file, _param.threads, std::nullopt);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  // Filter the input VCF in parallel
  ParallelFiltering(out_file);

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::FilterTumorOnlyTe() {
  VerifyModelCompatibility(_param.model, _model_config.scoring_cols);

  // load the VCF File
  Logging::Info("Loading VCF file {}", _param.vcf_file);
  io::VcfReader reader(_param.vcf_file);
  if (1 != reader.GetHeader()->GetNumSamples()) {
    throw error::Error("VCF file must contain exactly one sample");
  }

  CreateParentDirectoryIfNotExists(_param.vcf_output);

  Logging::Info("Writing output to VCF file {}", _param.vcf_output.string());

  _hdr = reader.GetHeader()->Clone();
  if (_param.command_line) {
    _hdr->AddCustomMetaDataLine(GetVcfHeaderLine(*_param.command_line));
  }
  // Append descriptions for the FILTER fields
  _hdr->AddFilterLine({kFilteringPassId, kFilteringPassDesc});
  _hdr->AddFilterLine({kFilteringMapQualityId, kFilteringMapQualityDesc});
  _hdr->AddFilterLine({kFilteringBlocklistedId, kFilteringBlocklistedDesc});
  _hdr->AddFilterLine({kFilteringMLScoreId, kFilteringMLScoreDesc});
  _hdr->AddFilterLine({kFilteringCountsId, kFilteringCountsDesc});
  _hdr->AddFilterLine({kFilteringBaseQualityId, kFilteringBaseQualityDesc});
  _hdr->AddFilterLine({kFilteringForcedId, kFilteringForcedDesc});
  _hdr->AddFilterLine({kFilteringAFId, kFilteringAFDesc});
  _hdr->AddFilterLine({kFilteringMinAltCountsId, kFilteringMinAltCountsDesc});

  // Append descriptions for the FORMAT fields
  _hdr->AddFormatLine({kWeightedCountsId, kWeightedCountsDesc, io::kNumberEachAllele, io::FieldType::kFloat});
  _hdr->AddFormatLine({kMachineLearningId, kMachineLearningDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kNonDuplexCountsId, kNonDuplexCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kDuplexCountsId, kDuplexCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kStrandBiasMetricId, kStrandBiasMetricDesc, io::kNumberOne, io::FieldType::kFloat});
  _hdr->AddFormatLine({kPlusOnlyCountsId, kPlusOnlyCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kMinusOnlyCountsId, kMinusOnlyCountsDesc, io::kNumberOne, io::FieldType::kInteger});
  _hdr->AddFormatLine({kSequenceContextId, kSequenceContextDesc, io::kNumberOne, io::FieldType::kString});
  _hdr->AddFormatLine({kPhysicalPhasedId, kPhysicalPhasedDesc, io::kNumberOne, io::FieldType::kString});
  _hdr->AddFormatLine({kHotspotId, kHotspotDesc, io::kNumberOne, io::FieldType::kInteger});

  const io::VcfWriter out_file(_param.vcf_output, _hdr);
  out_file.WriteHeader();

  // Load BED regions if provided
  const auto bed_regions = GetChromIntervalMap(_param.bed_file);
  // Partition the VCF file into regions for parallel processing
  _partitioned_regions = PartitionVcfRegions(_param.vcf_file, _param.threads, bed_regions);

  // Load reference sequences based on which chromosomes have variants in the VCF that need to be processed
  SetReferenceSequences();

  // Filter E2E by region
  ParallelFiltering(out_file);

  Logging::Info("Wrote output VCF file {}", _param.vcf_output);
}

void FilterVariantsClass::CheckTumorNormalSampleName() const {
  if (!_param.tumor_sample_name.has_value() || _param.tumor_sample_name->empty()) {
    throw error::Error("Tumor sample name for read groups is not specified");
  }
  Logging::Info("Tumor sample name for read groups: {}", _param.tumor_sample_name.value());
}

void FilterVariantsClass::VerifyParameters() const {
  using enum Workflow;
  switch (_model_config.workflow) {
    case kGermlineTagging:
    case kTumorNormalWgs: {
      CheckTumorNormalSampleName();
      break;
    }
    case kTumorOnlyTe:
    case kGermline:
    case kGermlineMultiSample:
      break;
    default:
      throw error::Error("Unsupported workflow");
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
  using enum Workflow;
  // This function acts as the entry point for the filter-variants tool. Based on the workflow the respective filtering
  // process is run

  // Create Filtering object
  FilterVariantsClass filtering(param);

  // Verify that the correct number of models have been passed based on the workflow, and that BAM input has been passed
  filtering.VerifyParameters();

  switch (param.workflow) {
    case kGermlineMultiSample:
    case kGermline: {
      filtering.FilterGermline();
      break;
    }
    case kGermlineTagging: {
      filtering.FilterGermlineTagging();
      break;
    }
    case kTumorOnlyTe: {
      filtering.FilterTumorOnlyTe();
      break;
    }
    case kTumorNormalWgs: {
      filtering.FilterTumorNormal();
      break;
    }
    default:
      throw error::Error("Specified Workflow has not yet been implemented");
  }
  Logging::Info("Filtered all variants");
}

}  // namespace xoos::svc
