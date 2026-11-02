#include "cli/define-copy-number-caller-options.h"

#include <cstddef>
#include <string>

#include <fmt/format.h>

#include <xoos/cli/cli.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/bed-validator.h>
#include <xoos/cli/validators/file-permission-validator.h>
#include <xoos/cli/validators/indexed-bam-validator.h>
#include <xoos/cli/validators/indexed-fasta-validator.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>

#include "CLI/CLI.hpp"
#include "cli/copy-number-caller-cli-option-names.h"
#include "cli/copy-number-caller-options.h"
#include "cli/vcf-parsing-options.h"
#include "coverage/calculate-coverage-options.h"
#include "intervals/augment-intervals-options.h"
#include "likelihood/likelihood-options.h"
#include "purity-ploidy-search/purity-ploidy-search-options.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {

namespace opt = cli_opt_name;

void DefineSomaticCommonOptions(CLI::App* const app,
                                CopyNumberCallerOptions& options,
                                const SegmentUndoMode undo_mode) {
  const bool hide = true;
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options.output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options.overwrite);
  AddOutputPrefixOption(app, opt::kOutputOptionsGroup, options.output_prefix);
  DefineIntervalOptions(app,
                        options.augment_intervals_options,
                        hide,
                        kAugmentIntervalsDefaultWholeGenomeIntervalSize,
                        kAugmentIntervalsDefaultWholeGenomeMinMappability);
  AddSexSelectionOption(app, opt::kSampleMetadataOptionsGroup, options.sample_metadata_options.sex_selection);
  AddTumorSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options.sample_metadata_options.tumor_sample_name);
  app->add_option(opt::kTumorBam,
                  options.tumor_bam_fname,
                  "Path to the input tumor BAM file. The BAM file must be coordinate-sorted and indexed.")
      ->check(cli::IndexedBamFileValidator())
      ->required()
      ->group(opt::kInputOptionsGroup);
  DefineAugmentIntervalsInputOptions(app, options, !hide, false);
  AddSeedSegmentsOption(app, opt::kInputOptionsGroup, options.seed_segments_fname);
  DefineCalculateCoverageAlgorithmOptions(app, options.calculate_coverage_options, hide);
  DefineGCCorrectAlgorithmOptions(app, options.gc_correct_options, hide);
  DefineSegmentationAlgorithmOptions(app, options.segmentation_options, hide, undo_mode);
  // BAF distance weight only applies when a BAF channel is present (somatic tumor/normal), so it is
  // registered here rather than in the shared DefineSegmentationAlgorithmOptions.
  AddBafDistanceWeightOption(app, opt::kSegmentationOptionsGroup, options.segmentation_options.baf_distance_weight);
  DefineSomaticSegmentationOptions(app, options.segmentation_options, hide);
  DefinePurityPloidySearchAlgorithmOptions(app, options.purity_ploidy_search_options, hide);
  DefineCoverageOptions(app, options);
  DefineBafFilteringOptions(app, options, hide);
  AddSaveTaskflowGraphOption(app, opt::kHiddenOptionsGroup, options.save_taskflow_graph);
  cli::AddThreadCountOption(app, opt::kThreads, options.threads)->group(opt::kPerformanceOptionsGroup);
}

void DefineSomaticTumorNormalWGSOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  app->add_option(opt::kNormalBam,
                  options->normal_bam_fname,
                  "Path to the input normal BAM file. The BAM file must be coordinate-sorted and indexed.")
      ->check(cli::IndexedBamFileValidator())
      ->required()
      ->group(opt::kInputOptionsGroup);
  // Register --vcf before DefineSomaticCommonOptions so that VCF-dependent options (e.g. baf filter
  // options that ->needs("--vcf")) can find it. The required() and needs() constraints are applied
  // after --normal-sample-name and --tumor-sample-name are registered.
  auto* vcf_opt = AddVcfOption(app,
                               opt::kInputOptionsGroup,
                               options->vcf_fname,
                               "Path to an input VCF file with germline variant sites. Heterozygous germline calls "
                               "allow for estimation of allele-specific copy number.")
                      ->required();
  DefineSomaticCommonOptions(app, *options, SegmentUndoMode::kIterative);
  // Register --normal-sample-name after DefineSomaticCommonOptions (which registers --tumor-sample-name)
  // so vcf_opt can declare needs on both sample name options.
  AddNormalSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.normal_sample_name);
  vcf_opt->needs(opt::kNormalSampleName)->needs(opt::kTumorSampleName);
  // We also include augment-intervals options here because tumor-normal will require an
  // augmented-intervals file created from the whole genome. We do not require the user to run AugmentIntervals
  // themselves for tumor-normal. (for TE, AugmentIntervals is run as a part of GeneratePanelOfNormals)
  AddIntervalsOption(app, opt::kHiddenOptionsGroup, options->augmented_intervals_fname);
  AddSaveBafSegmentsOption(app, opt::kHiddenOptionsGroup, options->save_baf_segments);
  // Register --tumor-ploidy before --tumor-purity, then wire up mutual needs after both are registered
  auto* ploidy_opt =
      AddTumorPloidyOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.ploidy);
  auto* purity_opt =
      AddTumorPurityOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.purity);
  purity_opt->needs(ploidy_opt);
  ploidy_opt->needs(purity_opt);
}

void DefineSomaticTumorOnlyWGSOptions(CLI::App* const app, const CopyNumberCallerOptionsPtr& options) {
  // hide advanced algorithm options from the default --help listing
  const bool hide = true;

  // --- Input: tumor BAM (no matched normal for tumor-only) ---
  app->add_option(opt::kTumorBam,
                  options->tumor_bam_fname,
                  "Path to the input tumor BAM file. The BAM file must be coordinate-sorted and indexed.")
      ->check(cli::IndexedBamFileValidator())
      ->required()
      ->group(opt::kInputOptionsGroup);

  // --- Output ---
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  AddOutputPrefixOption(app, opt::kOutputOptionsGroup, options->output_prefix);

  // --- Intervals / augment-intervals inputs (whole genome, generated internally) ---
  DefineIntervalOptions(app,
                        options->augment_intervals_options,
                        hide,
                        kAugmentIntervalsDefaultWholeGenomeIntervalSize,
                        kAugmentIntervalsDefaultWholeGenomeMinMappability);
  DefineAugmentIntervalsInputOptions(app, *options, !hide, false);
  AddSeedSegmentsOption(app, opt::kInputOptionsGroup, options->seed_segments_fname);
  AddIntervalsOption(app, opt::kHiddenOptionsGroup, options->augmented_intervals_fname);

  // --- Sample metadata ---
  AddSexSelectionOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex_selection);
  AddTumorSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.tumor_sample_name)
      ->required();

  // --- Tumor coverage / GC correction ---
  AddTumorMinCoverageOption(app, opt::kReadFilteringOptionsGroup, options->tumor_min_coverage);
  DefineCalculateCoverageAlgorithmOptions(app, options->calculate_coverage_options, hide);
  DefineGCCorrectAlgorithmOptions(app, options->gc_correct_options, hide);

  // --- Segmentation (somatic, iterative log-ratio undoing) ---
  DefineSegmentationAlgorithmOptions(app, options->segmentation_options, hide, SegmentUndoMode::kIterative);
  AddBafDistanceWeightOption(app, opt::kSegmentationOptionsGroup, options->segmentation_options.baf_distance_weight);
  DefineSomaticSegmentationOptions(app, options->segmentation_options, hide);

  // --- Purity / ploidy search ---
  DefinePurityPloidySearchAlgorithmOptions(
      app, options->purity_ploidy_search_options, hide, kTumorOnlyWgsBafLikelihoodWeight);
  app->get_option(opt::kBafLikelihoodWeightForPurityPloidySearchRerun)->group(opt::kPurityPloidySearchOptionsGroup);

  // --- Tumor-only BAF extraction (population VCF drives germline het identification) ---
  AddVcfOption(app,
               opt::kInputOptionsGroup,
               options->vcf_fname,
               "Path to the single-sample tumor VCF file with variant sites. Must be position-sorted.")
      ->required();
  AddVcfOption(app,
               opt::kInputOptionsGroup,
               options->population_vcf_fname,
               "Path to the population SNP VCF (e.g., gnomAD). Used to identify likely germline heterozygous sites. "
               "Must be indexed with .tbi.",
               opt::kPopulationVcf)
      ->required();
  app->add_option(opt::kSnpBlocklistBed,
                  options->snp_blocklist_bed_fname,
                  "Path to a BED file of 0-based regions whose overlapping SNP sites are excluded from BAF "
                  "extraction (e.g., known problematic or recurrent-artifact loci). When omitted, uses the bundled "
                  "hg38 SNP blocklist at /resources/hg38_snp_blocklist.bed.")
      ->check(cli::BedFileValidator())
      ->group(opt::kInputOptionsGroup);
  app->add_option(opt::kMinPopulationAf,
                  options->min_population_af,
                  "Minimum population allele frequency to consider a site as a likely germline het.")
      ->default_val(kVcfParsingOptionsDefaultMinPopulationAF)
      ->check(CLI::Range(0.0, 1.0))
      ->group(opt::kVcfParsingOptionsGroup);
  app->add_option(
         opt::kMinTumorDepth, options->min_tumor_depth, "Minimum total depth (ref + alt) at a site in the tumor VCF.")
      ->default_val(kVcfParsingOptionsDefaultTumorOnlyMinDepth)
      ->check(CLI::PositiveNumber)
      ->group(opt::kVcfParsingOptionsGroup);
  app->add_option(
         opt::kMinTumorBaf, options->baf_filter_options.tumor_sample_min_baf, "Minimum BAF to include in the output.")
      ->default_val(kVcfParsingOptionsDefaultTumorOnlyMinBAF)
      ->check(CLI::Range(0.0, 1.0))
      ->group(opt::kVcfParsingOptionsGroup);
  app->add_option(
         opt::kMaxTumorBaf, options->baf_filter_options.tumor_sample_max_baf, "Maximum BAF to include in the output.")
      ->default_val(kVcfParsingOptionsDefaultTumorOnlyMaxBAF)
      ->check(CLI::Range(0.0, 1.0))
      ->group(opt::kVcfParsingOptionsGroup);
  AddTumorSampleMinDepthOption(app, opt::kHiddenOptionsGroup, options->baf_filter_options.tumor_sample_min_depth);
  AddAdInOption(app, opt::kHiddenOptionsGroup, options->ad_fname);

  // --- Diagnostics / performance ---
  AddSaveBafSegmentsOption(app, opt::kHiddenOptionsGroup, options->save_baf_segments);
  AddSaveTaskflowGraphOption(app, opt::kHiddenOptionsGroup, options->save_taskflow_graph);
  cli::AddThreadCountOption(app, opt::kThreads, options->threads)->group(opt::kPerformanceOptionsGroup);

  // --- Optional user-provided purity/ploidy (mutually dependent) ---
  auto* const ploidy_opt =
      AddTumorPloidyOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.ploidy);
  auto* const purity_opt =
      AddTumorPurityOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.purity);
  purity_opt->needs(ploidy_opt);
  ploidy_opt->needs(purity_opt);
}

void DefineGermlineNormalWGSOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  // hide options for advanced users
  const bool hide = true;
  app->add_option(opt::kBam,
                  options->normal_bam_fname,
                  "Path to the input BAM file. The BAM file must be coordinate-sorted and indexed.")
      ->check(cli::IndexedBamFileValidator())
      ->required()
      ->group(opt::kInputOptionsGroup);
  AddLibraryTypeOption(app, opt::kInputOptionsGroup, options->library_type);
  DefineAugmentIntervalsInputOptions(app, *options, false, false);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  AddOutputPrefixOption(app, opt::kOutputOptionsGroup, options->output_prefix);
  DefineIntervalOptions(
      app,
      options->augment_intervals_options,
      hide,
      kAugmentIntervalsGermlineWGSDefaultIntervalSize,
      kAugmentIntervalsGermlineWGSDefaultWholeGenomeMinMappability,
      fmt::format("If not specified, defaults to {} bp for SBX-S inputs and {} bp otherwise, based on the resolved "
                  "`--library-type`.",
                  kAugmentIntervalsGermlineWGSSbxSIntervalSize,
                  kAugmentIntervalsGermlineWGSDefaultIntervalSize));
  AddSeedSegmentsOption(app, opt::kInputOptionsGroup, options->seed_segments_fname);
  // Register --sample-name before --vcf so VCF can declare a dependency on it
  AddSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.normal_sample_name);
  AddVcfOption(app,
               opt::kInputOptionsGroup,
               options->vcf_fname,
               "Path to an input VCF file for germline SNPs (for BAF visualization).")
      ->needs(opt::kSampleName);
  DefineCalculateCoverageAlgorithmOptions(app, options->calculate_coverage_options, hide);
  DefineSegmentationAlgorithmOptions(app, options->segmentation_options, hide, SegmentUndoMode::kDisabled);
  DefineGermlineCallFlagOptions(app, options->likelihood_options, hide);
  AddSexSelectionOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex_selection);
  cli::AddThreadCountOption(app, opt::kThreads, options->threads)->group(opt::kPerformanceOptionsGroup);
  // Hidden advanced options
  AddNormalSampleMinDepthOption(app,
                                opt::kHiddenOptionsGroup,
                                options->baf_filter_options.normal_sample_min_depth,
                                kVcfParsingOptionsDefaultGermlineNormalSampleMinDepth);
  AddGcCorrectFirstSpanOption(app, opt::kHiddenOptionsGroup, options->gc_correct_options.first_span);
  AddIntervalsOption(app, opt::kHiddenOptionsGroup, options->augmented_intervals_fname);
  AddSaveLogRatioSegmentsOption(app, opt::kHiddenOptionsGroup, options->save_log_ratio_segments);
  AddSaveTaskflowGraphOption(app, opt::kHiddenOptionsGroup, options->save_taskflow_graph);
}

// ==============================================================================
// AugmentIntervals
// ==============================================================================
void DefineAugmentIntervalsInputOptions(CLI::App* const app,
                                        CopyNumberCallerOptions& options,
                                        const bool hide,
                                        const bool require_mappability) {
  app->add_option(opt::kReference,
                  options.reference_genome_fname,
                  "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->check(cli::IndexedFastaFileValidator())
      ->group(opt::kInputOptionsGroup)
      ->required();
  auto* const mappability_option =
      AddMappabilityBigwigOption(app, opt::kInputOptionsGroup, options.mappability_bigwig_fname);
  if (require_mappability) {
    mappability_option->required();
  } else {
    mappability_option->description(
        "Path to an input BigWig file describing mappability across the genome. When omitted, SBX-S input uses "
        "/resources/hg38_500mer_mappability.bw and other or unresolved library types use "
        "/resources/hg38_mappability.bw.");
  }
  AddBlocklistBedOption(app, hide ? opt::kHiddenOptionsGroup : opt::kInputOptionsGroup, options.blocklist_bed_fname);
}

void DefineIntervalOptions(CLI::App* const app,
                           AugmentIntervalsOptions& options,
                           const bool hide,
                           const size_t default_interval_size,
                           const f64 default_mappability,
                           const std::string& interval_size_description_suffix) {
  AddWholeGenomeIntervalSizeOption(app,
                                   opt::kRegionOptionsGroup,
                                   options.whole_genome_interval_size,
                                   default_interval_size,
                                   interval_size_description_suffix);
  AddWholeGenomeMinMappabilityOption(
      app, opt::kRegionOptionsGroup, options.whole_genome_min_mappability, default_mappability);
  AddMinGcContentOption(app, hide ? opt::kHiddenOptionsGroup : opt::kRegionOptionsGroup, options.min_gc_content);
  AddMaxGcContentOption(app, hide ? opt::kHiddenOptionsGroup : opt::kRegionOptionsGroup, options.max_gc_content);
}

void DefineIntervalOptions(CLI::App* const app,
                           AugmentIntervalsOptions& options,
                           const bool hide,
                           const size_t default_interval_size,
                           const f64 default_mappability) {
  DefineIntervalOptions(app, options, hide, default_interval_size, default_mappability, "");
}

void DefineAdvancedIntervalOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    AugmentIntervalsOptions& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddNoOffTargetsOption(app, opt::kRegionOptionsGroup, options.no_off_targets);
  AddOnTargetIntervalSizeOption(app, opt::kRegionOptionsGroup, options.on_target_interval_size);
  AddOffTargetMinWidthOption(app, opt::kRegionOptionsGroup, options.off_target_min_width);
  AddOffTargetTrimLengthOption(app, opt::kRegionOptionsGroup, options.off_target_trim_length);
  AddOffTargetIntervalSizeOption(app, opt::kRegionOptionsGroup, options.off_target_interval_size);
  AddOnTargetMinMappabilityOption(app, opt::kRegionOptionsGroup, options.on_target_min_mappability);
  AddOffTargetMinMappabilityOption(app, opt::kRegionOptionsGroup, options.off_target_min_mappability);
  AddSexChromosomeTargetMinMappabilityOption(app, opt::kRegionOptionsGroup, options.sex_chromosome_min_mappability);
}

void DefineAugmentIntervalsSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  DefineAugmentIntervalsInputOptions(app, *options, false, true);
  app->add_option(opt::kBam,
                  options->normal_bam_fname,
                  "Path to the input BAM file. Used to check whether the chrY PAR region was masked during "
                  "alignment. If unmasked, PAR intervals are removed from the output intervals file.")
      ->check(cli::IndexedBamFileValidator())
      ->required()
      ->group(opt::kInputOptionsGroup);
  AddSeedSegmentsOption(app, opt::kInputOptionsGroup, options->seed_segments_fname)->required();
  AddBedOption(app, opt::kInputOptionsGroup, options->augmented_intervals_fname);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  AddWholeGenomeOption(
      app, opt::kRegionOptionsGroup, options->augment_intervals_options.generate_whole_genome_intervals);
  DefineIntervalOptions(app,
                        options->augment_intervals_options,
                        false,
                        kAugmentIntervalsDefaultWholeGenomeIntervalSize,
                        kAugmentIntervalsDefaultWholeGenomeMinMappability);
  DefineAdvancedIntervalOptions(app, options->augment_intervals_options);
  AddIncludeAltContigsOption(app, opt::kRegionOptionsGroup, options->augment_intervals_options.include_alt_contigs);
}

// ==============================================================================
// CalculateCoverage
// ==============================================================================
void DefineCalculateCoverageAlgorithmOptions(CLI::App* const app, CalculateCoverageOptions& options, const bool hide) {
  // we are not validating this as an s32 yet; currently validation happens with an htslib call down stream
  AddCoverageExcludeFlagsOption(app, opt::kReadFilteringOptionsGroup, options.exclude_flags);
  AddCoverageIgnoreDnOption(
      app, hide ? opt::kHiddenOptionsGroup : opt::kCalculateCoverageOptionsGroup, options.ignore_DN);
}

void DefineCalculateCoverageSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  app->add_option(opt::kNormalBam,
                  options->normal_bam_fname,
                  "Path to the input BAM file. The BAM file must be coordinate-sorted and indexed.")
      ->check(cli::IndexedBamFileValidator())
      ->required()
      ->group(opt::kInputOptionsGroup);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  AddIntervalsOption(app, opt::kInputOptionsGroup, options->augmented_intervals_fname)->required();
  DefineCalculateCoverageAlgorithmOptions(app, options->calculate_coverage_options, false);
  cli::AddThreadCountOption(app, opt::kThreads, options->threads)->group(opt::kPerformanceOptionsGroup);
}

// ==============================================================================
// GCCorrect
// ==============================================================================
void DefineGCCorrectAlgorithmOptions(CLI::App* const app, GCCorrectOptions& options, const bool hide) {
  const auto* const group_name = hide ? opt::kHiddenOptionsGroup : opt::kGCCorrectOptionsGroup;
  AddGcCorrectFirstSpanOption(app, group_name, options.first_span);
}

void DefineGCCorrectSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddCoverageFileOption(app, opt::kInputOptionsGroup, options->normal_coverage_fname)->required();
  AddIntervalsOption(app, opt::kInputOptionsGroup, options->augmented_intervals_fname)->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  DefineGCCorrectAlgorithmOptions(app, options->gc_correct_options, false);
  cli::AddThreadCountOption(app, opt::kThreads, options->threads)->group(opt::kPerformanceOptionsGroup);
}

// ==============================================================================
// TwoSampleLogR
// ==============================================================================
void DefineCoverageOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptions& options) {  // NOSONAR cpp:S995 cpp:M23_131
                                         // by CLI11 on parse; non-const required
  AddTumorMinCoverageOption(app, opt::kReadFilteringOptionsGroup, options.tumor_min_coverage);
  AddNormalMinCoverageOption(app, opt::kReadFilteringOptionsGroup, options.normal_min_coverage);
}

void DefineTwoSampleLogRSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddTumorCoverageFileOption(app, opt::kInputOptionsGroup, options->tumor_coverage_fname)->required();
  AddNormalCoverageFileOption(app, opt::kInputOptionsGroup, options->normal_coverage_fname)->required();
  app->add_option(opt::kReference,
                  options->reference_genome_fname,
                  "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->check(cli::IndexedFastaFileValidator())
      ->group(opt::kInputOptionsGroup);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  DefineCoverageOptions(app, *options);
}

// ==============================================================================
// OneSampleLogR
// ==============================================================================
void DefineOneSampleLogRSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddCoverageFileOption(app, opt::kInputOptionsGroup, options->normal_coverage_fname)->required();
  app->add_option(opt::kReference,
                  options->reference_genome_fname,
                  "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->check(cli::IndexedFastaFileValidator())
      ->group(opt::kInputOptionsGroup);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
}

// ==============================================================================
// Segmentation
// ==============================================================================
void DefineSegmentationAlgorithmOptions(CLI::App* const app,
                                        SegmentationOptions& options,
                                        const bool hide,
                                        const SegmentUndoMode undo_mode) {
  const auto* const group_name = hide ? opt::kHiddenOptionsGroup : opt::kSegmentationOptionsGroup;
  AddPruningClusteringParameterOption(app, opt::kSegmentationOptionsGroup, options.pruning_clustering_parameter);
  AddUndoLogrSegmentsOption(app, opt::kSegmentationOptionsGroup, options.undo_logr_mode, undo_mode);
  AddUndoLogrSegmentsSdFactorOption(app, opt::kSegmentationOptionsGroup, options.undo_logr_segments_sd_factor);
  AddMaxNumSegmentsOption(app, opt::kSegmentationOptionsGroup, options.max_num_segments);
  AddIncrementUndoLogrSegmentsSdFactorOption(app, group_name, options.increment_undo_logr_segments_sd_factor);
  AddDisableMergingOption(app, group_name, options.disable_merging);
  AddDisableHierarchicalPruningOption(app, group_name, options.disable_hierarchical_pruning);
}

void DefineSomaticSegmentationOptions(CLI::App* const app, SegmentationOptions& options, const bool hide) {
  const auto* const group_name = hide ? opt::kHiddenOptionsGroup : opt::kSegmentationOptionsGroup;
  AddMinTForAutomaticSegmentationOption(app, group_name, options.min_t_for_automatic_segmentation);
  AddUndoDhSegmentsOption(app, group_name, options.undo_dh_segments);
  AddUndoDhSegmentsSdFactorOption(app, group_name, options.undo_dh_segments_sd_factor)
      ->needs(opt::kEnableDhSegmentUndoing);
}

void DefineBafFilteringOptions(CLI::App* const app, CopyNumberCallerOptions& options, const bool hide) {
  const auto* const group_name = hide ? opt::kHiddenOptionsGroup : opt::kVcfParsingOptionsGroup;
  AddNormalSampleMinDepthOption(app,
                                opt::kVcfParsingOptionsGroup,
                                options.baf_filter_options.normal_sample_min_depth,
                                kVcfParsingOptionsDefaultSomaticNormalSampleMinDepth);
  AddNormalSampleMinBafOption(app, opt::kVcfParsingOptionsGroup, options.baf_filter_options.normal_sample_min_baf);
  AddNormalSampleMaxBafOption(app, opt::kVcfParsingOptionsGroup, options.baf_filter_options.normal_sample_max_baf);
  AddTumorSampleMinDepthOption(app, group_name, options.baf_filter_options.tumor_sample_min_depth);
  AddTumorSampleMinBafOption(app, group_name, options.baf_filter_options.tumor_sample_min_baf);
  AddTumorSampleMaxBafOption(app, group_name, options.baf_filter_options.tumor_sample_max_baf);
  AddAdInOption(app, group_name, options.ad_fname);
  AddForceEnableSomaticVariantParsingOption(
      app, opt::kVcfParsingOptionsGroup, options.baf_filter_options.force_enable_somatic_variant_parsing);
}

void DefineSegmentationSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddLogrsOption(app, opt::kInputOptionsGroup, options->logrs_fname)->required();
  AddSeedSegmentsOption(app, opt::kInputOptionsGroup, options->seed_segments_fname)->required();
  // Register --vcf before --bafs and --reference so exclusion constraints can reference it
  auto* vcf_opt = AddVcfOption(
      app, opt::kInputOptionsGroup, options->vcf_fname, "Path to an input VCF file with germline variant sites.");
  // --bafs: pre-computed BAF BED file, mutually exclusive with --vcf and --ad
  auto* bafs_opt = AddBafsOption(app, opt::kInputOptionsGroup, options->bafs_fname);
  bafs_opt->excludes(vcf_opt);
  // --reference is optional; allowed with either --vcf or --bafs (for BigWig output)
  app->add_option(opt::kReference,
                  options->reference_genome_fname,
                  "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->check(cli::IndexedFastaFileValidator())
      ->group(opt::kInputOptionsGroup);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  AddSegmentationModeOption(app, opt::kSegmentationMode, options->segmentation_options.segmentation_mode)->required();
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddNormalSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.normal_sample_name);
  AddTumorSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.tumor_sample_name);
  AddMinObservationsPerSegmentOption(
      app, opt::kSegmentationOptionsGroup, options->segmentation_options.min_obs_per_segment);
  DefineSegmentationAlgorithmOptions(app, options->segmentation_options, false, SegmentUndoMode::kIterative);
  // BAF distance weight only applies when a BAF channel is present, so it is registered here rather
  // than in the shared DefineSegmentationAlgorithmOptions (which is also used by germline-normal-wgs).
  AddBafDistanceWeightOption(app, opt::kSegmentationOptionsGroup, options->segmentation_options.baf_distance_weight);
  DefineSomaticSegmentationOptions(app, options->segmentation_options, false);
  AddCbsMethodOption(app, opt::kSegmentationOptionsGroup, options->segmentation_options.cbs_method);
  AddTruncateOutliersOption(app, opt::kSegmentationOptionsGroup, options->segmentation_options.truncate_outliers);
  // Register the VCF-parsing (module-specific) options before the Performance option so the help menu
  // renders module-specific sections before "Performance Options" (CLI11 orders groups by first use).
  DefineBafFilteringOptions(app, *options, false);
  cli::AddThreadCountOption(app, opt::kThreads, options->threads)->group(opt::kPerformanceOptionsGroup);

  // Wire mutual exclusivity: --bafs excludes --ad and all VCF-parsing options
  for (const auto* excluded : {opt::kAlleleDepthTsv,
                               opt::kMinVcfNormalDepth,
                               opt::kMinNormalBaf,
                               opt::kMaxNormalBaf,
                               opt::kMinVcfTumorDepth,
                               opt::kMinTumorBaf,
                               opt::kMaxTumorBaf,
                               opt::kForceEnableSomaticVariantParsing}) {
    app->get_option(excluded)->excludes(bafs_opt);
  }
}

// ==============================================================================
// PurityPloidySearch
// ==============================================================================
void DefinePurityPloidySearchAlgorithmOptions(CLI::App* const app,
                                              PurityPloidySearchOptions& options,
                                              const bool hide) {
  DefinePurityPloidySearchAlgorithmOptions(app, options, hide, kBafLikelihoodWeight);
}

void DefinePurityPloidySearchAlgorithmOptions(CLI::App* const app,
                                              PurityPloidySearchOptions& options,
                                              const bool hide,
                                              const f64 baf_likelihood_weight_default) {
  const auto* const group_name = hide ? opt::kHiddenOptionsGroup : opt::kPurityPloidySearchOptionsGroup;
  AddSegMaxLogrOption(app, group_name, options.seg_max_logr);
  AddSegMinNumLogrsOption(app, group_name, options.seg_min_num_logrs);
  AddSegMinNumSnpsOption(app, group_name, options.seg_min_num_snps);
  AddBafLikelihoodWeightForPurityPloidySearchRerunOption(
      app, group_name, options.baf_likelihood_weight, baf_likelihood_weight_default);
  AddMinLohFractionForWgdOption(app, opt::kHiddenOptionsGroup, options.min_loh_fraction_for_wgd);
  // Diagnostic output flag is always hidden regardless of the algorithm-option grouping.
  AddWriteSegmentLikelihoodOption(app, opt::kHiddenOptionsGroup, options.write_segment_likelihood);
}

void DefinePurityPloidySearchSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddSegmentsOption(app, opt::kInputOptionsGroup, options->segments_fname)->required();
  AddLogrsOption(app, opt::kInputOptionsGroup, options->logrs_fname)->required();
  AddBafsOption(app, opt::kInputOptionsGroup, options->bafs_fname)->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddSampleIdOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sample_id);
  DefinePurityPloidySearchAlgorithmOptions(app, options->purity_ploidy_search_options, false);
  AddTumorOnlyOption(app,
                     opt::kPurityPloidySearchOptionsGroup,
                     options->purity_ploidy_search_options.tumor_only,
                     "Run the purity/ploidy search in tumor-only mode. Segments are classified as "
                     "germline or somatic so that rare germline CNVs do not bias the purity estimate, and allosomes "
                     "are excluded from the search.");
}

// ==============================================================================
// Likelihood (PredictSomaticCNA / PredictGermlineCNV)
// ==============================================================================
void DefinePredictSomaticCNASubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddSegmentsOption(app, opt::kInputOptionsGroup, options->segments_fname)->required();
  AddLogrsOption(app, opt::kInputOptionsGroup, options->logrs_fname)->required();
  AddBafsOption(app, opt::kInputOptionsGroup, options->bafs_fname)->required();
  app->add_option(opt::kReference,
                  options->reference_genome_fname,
                  "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->check(cli::IndexedFastaFileValidator())
      ->group(opt::kInputOptionsGroup)
      ->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  AddOutputPrefixOption(app, opt::kOutputOptionsGroup, options->output_prefix);
  // Register --tumor-ploidy before --tumor-purity, then wire up mutual needs after both are registered
  auto* ploidy_opt =
      AddTumorPloidyOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.ploidy)->required();
  auto* purity_opt =
      AddTumorPurityOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.purity)->required();
  purity_opt->needs(ploidy_opt);
  ploidy_opt->needs(purity_opt);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddTumorOnlyOption(app,
                     opt::kInputOptionsGroup,
                     options->likelihood_options.tumor_only,
                     "Predict copy number in tumor-only mode. Segments are classified as germline "
                     "or somatic; germline segments are flagged with the GERMLINE INFO field in the output VCF.");
}

void DefineGermlineCallFlagOptions(CLI::App* const app, LikelihoodOptions& options, const bool hide) {
  AddMapqCutoffToFlagCallsOption(app, opt::kOutputOptionsGroup, options.mapq_cutoff_for_calls);
  AddCnvLengthFlagMinSizeOption(
      app, hide ? opt::kHiddenOptionsGroup : opt::kOutputOptionsGroup, options.cnv_length_flag_min_size);
}

void DefinePredictGermlineCNVSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddSegmentsOption(app, opt::kInputOptionsGroup, options->segments_fname)->required();
  AddLogrsOption(app, opt::kInputOptionsGroup, options->logrs_fname)->required();
  AddMapqsOption(app, opt::kInputOptionsGroup, options->mapping_qualities_fname);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddSampleIdOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sample_id);
  DefineGermlineCallFlagOptions(app, options->likelihood_options, false);
}

// ==============================================================================
// MergeSegments
// ==============================================================================
void DefineMergeSegmentsSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddSegmentsOption(app, opt::kInputOptionsGroup, options->segments_fname)->required();
  AddLogrsOption(app, opt::kInputOptionsGroup, options->logrs_fname)->required();
  AddMapqsOption(app, opt::kInputOptionsGroup, options->mapping_qualities_fname)->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  AddUseMapqsFromFileOption(app, opt::kInputOptionsGroup, options->merge_segments_options.use_mapqs_observations);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddSampleIdOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sample_id);
  AddMinMapqOption(app, opt::kOutputOptionsGroup, options->merge_segments_options.min_mapq_threshold);
  AddRecalculatePerSegmentDataOption(
      app, opt::kOutputOptionsGroup, options->merge_segments_options.recalculate_per_segment_data);
}

// ==============================================================================
// SegToVcf
// ==============================================================================
void DefineSegToVcfSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddSegmentsOption(app, opt::kInputOptionsGroup, options->segments_fname)->required();
  app->add_option(opt::kReference,
                  options->reference_genome_fname,
                  "Path to the input reference FASTA file. The FASTA file must be indexed with a .fai file.")
      ->check(cli::IndexedFastaFileValidator())
      ->group(opt::kInputOptionsGroup);
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
  AddSexOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sex)->required();
  AddSampleIdOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.sample_id);
}

void DefinePredictSexSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  app->add_option("--coverage-file", options->predict_sex_coverage_fname, "Coverage file (pre or post GC correction).")
      ->check(cli::FileReadableValidator())
      ->group(opt::kInputOptionsGroup)
      ->required();
  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);
}

// ==============================================================================
void DefineTumorOnlyBafsSubcommandOptions(
    CLI::App* const app,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CopyNumberCallerOptionsPtr& options) {  // NOSONAR cpp:S995 cpp:M23_131
  AddVcfOption(app,
               opt::kInputOptionsGroup,
               options->vcf_fname,
               "Path to the single-sample tumor VCF file. Must be position-sorted.")
      ->required();

  AddVcfOption(app,
               opt::kInputOptionsGroup,
               options->population_vcf_fname,
               "Path to the population SNP VCF (e.g., gnomAD). Used to identify likely germline heterozygous sites. "
               "Must be indexed with .tbi.",
               opt::kPopulationVcf)
      ->required();

  app->add_option(opt::kSnpBlocklistBed,
                  options->snp_blocklist_bed_fname,
                  "Path to a BED file of 0-based regions whose overlapping SNP sites are excluded from BAF "
                  "extraction (e.g., known problematic or recurrent-artifact loci).")
      ->check(cli::BedFileValidator())
      ->group(opt::kInputOptionsGroup);

  AddOutputDirOption(app, opt::kOutputOptionsGroup, options->output_dir);
  AddOverwriteOption(app, opt::kOutputOptionsGroup, options->overwrite);

  app->add_option(opt::kMinPopulationAf,
                  options->min_population_af,
                  "Minimum population allele frequency to consider a site as a likely germline het.")
      ->default_val(kVcfParsingOptionsDefaultMinPopulationAF)
      ->check(CLI::Range(0.0, 1.0))
      ->group(opt::kVcfParsingOptionsGroup);

  app->add_option(
         opt::kMinTumorDepth, options->min_tumor_depth, "Minimum total depth (ref + alt) at a site in the tumor VCF.")
      ->default_val(kVcfParsingOptionsDefaultTumorOnlyMinDepth)
      ->check(CLI::PositiveNumber)
      ->group(opt::kVcfParsingOptionsGroup);

  app->add_option(
         opt::kMinTumorBaf, options->baf_filter_options.tumor_sample_min_baf, "Minimum BAF to include in the output.")
      ->default_val(kVcfParsingOptionsDefaultTumorOnlyMinBAF)
      ->check(CLI::Range(0.0, 1.0))
      ->group(opt::kVcfParsingOptionsGroup);

  app->add_option(
         opt::kMaxTumorBaf, options->baf_filter_options.tumor_sample_max_baf, "Maximum BAF to include in the output.")
      ->default_val(kVcfParsingOptionsDefaultTumorOnlyMaxBAF)
      ->check(CLI::Range(0.0, 1.0))
      ->group(opt::kVcfParsingOptionsGroup);

  AddTumorSampleNameOption(app, opt::kSampleMetadataOptionsGroup, options->sample_metadata_options.tumor_sample_name)
      ->required();
}

}  // namespace xoos::cnc
