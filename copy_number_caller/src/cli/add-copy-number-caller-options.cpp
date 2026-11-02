#include <cstddef>
#include <optional>
#include <string>

#include <fmt/core.h>
#include <fmt/format.h>

#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/validators/bed-validator.h>
#include <xoos/cli/validators/file-permission-validator.h>
#include <xoos/cli/validators/filename-component-validator.h>
#include <xoos/cli/validators/indexed-vcf-validator.h>
#include <xoos/io/library-type.h>
#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "CLI/CLI.hpp"
#include "cli/copy-number-caller-cli-option-names.h"
#include "cli/copy-number-caller-options.h"
#include "cli/define-copy-number-caller-options.h"
#include "cli/vcf-parsing-options.h"
#include "coverage/calculate-coverage-options.h"
#include "gc-correct/gc-correct-options.h"
#include "intervals/augment-intervals-options.h"
#include "likelihood/likelihood-options.h"
#include "merge-segments/merge-segments-options.h"
#include "purity-ploidy-search/purity-ploidy-search-options.h"
#include "segmentation/max-t.h"
#include "segmentation/segmentation-mode.h"
#include "segmentation/segmentation-options.h"

namespace xoos::cnc {

namespace opt = cli_opt_name;

// =============================================================================
// Add*Option helper implementations
// =============================================================================
// Shared End-to-End options
CLI::Option* AddIntervalsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<fs::path>& binding) {  // NOSONAR cpp:S995 cpp:M23_131

  return app->add_option(opt::kIntervals, binding, "Annotated intervals file.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

// Shared IO options
CLI::Option* AddSegmentsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<fs::path>& binding) {  // NOSONAR cpp:S995 cpp:M23_131

  return app
      ->add_option(
          opt::kSegments, binding, "Initial copy number calls in SEG format, as output by the copy_number_caller.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddLogrsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<fs::path>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kLogRatios,
                   binding,
                   "Log ratios file in BED format (Contig, Start, End, LogR), as output by the copy_number_caller.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddBafsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<fs::path>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kBafs,
                   binding,
                   "Tumor BAF values in BED format (Contig, Start, End, BAF, Ref_AD, Alt_AD), as output by the "
                   "copy_number_caller. Used for estimation of allele-specific copy number.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddMapqsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<fs::path>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMapqs,
                   binding,
                   "Mapping quality file in BED format (Contig, Start, End, MeanMapQ). Use with "
                   "--mapq-cutoff-to-flag-calls.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddVcfOption(CLI::App* const app,
                          const std::string& group_name,
                          std::optional<fs::path>& binding,
                          const std::string& description) {
  // CLI11 on parse; non-const required
  return AddVcfOption(app, group_name, binding, description, opt::kVcf);
}

CLI::Option* AddVcfOption(CLI::App* const app,
                          const std::string& group_name,
                          std::optional<fs::path>& binding,
                          const std::string& description,
                          const std::string& option_name) {
  // CLI11 on parse; non-const required
  return app->add_option(option_name, binding, description)->check(cli::IndexedVcfFileValidator())->group(group_name);
}

// VCF Parsing
CLI::Option* AddNormalSampleMinDepthOption(CLI::App* const app,
                                           const std::string& group_name,
                                           s32& binding,
                                           const s32 default_depth) {
  return app
      ->add_option(opt::kMinVcfNormalDepth,
                   binding,
                   "Min normal depth for extracting variants from a matched tumor-normal or normal-only VCF. "
                   "Variants with depth less than this value are excluded.")
      ->check(CLI::TypeValidator<s32>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(default_depth)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddNormalSampleMinBafOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMinNormalBaf,
                   binding,
                   "Minimum BAF for a variant in the normal sample to be considered heterozygous. Variants below "
                   "this value are excluded.")
      ->check(CLI::TypeValidator<float>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(kVcfParsingOptionsDefaultSomaticNormalSampleMinBAF)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddNormalSampleMaxBafOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMaxNormalBaf,
                   binding,
                   "Maximum BAF for a variant in the normal sample to be considered heterozygous. Variants above "
                   "this value are excluded.")
      ->check(CLI::TypeValidator<float>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(kVcfParsingOptionsDefaultSomaticNormalSampleMaxBAF)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddTumorSampleMinDepthOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    s32& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kMinVcfTumorDepth, binding, "Minimum tumor depth for extracting variants.")
      ->check(CLI::TypeValidator<s32>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(kVcfParsingOptionsDefaultSomaticTumorSampleMinDepth)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddTumorSampleMinBafOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMinTumorBaf,
                   binding,
                   "Minimum BAF for a variant in the tumor sample to be extracted. Variants below this value are "
                   "excluded.")
      ->check(CLI::TypeValidator<float>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(kVcfParsingOptionsDefaultSomaticTumorSampleMinBAF)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddTumorSampleMaxBafOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMaxTumorBaf,
                   binding,
                   "Maximum BAF for a variant in the tumor sample to be extracted. Variants above this value are "
                   "excluded.")
      ->check(CLI::TypeValidator<float>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(kVcfParsingOptionsDefaultSomaticTumorSampleMaxBAF)
      ->needs(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddAdInOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<fs::path>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(
          opt::kAlleleDepthTsv, binding, "Alternative to VCF. Requires Chromosome, Position, and TumorBaf columns.")
      ->excludes(opt::kVcf)
      ->group(group_name);
}

CLI::Option* AddForceEnableSomaticVariantParsingOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kForceEnableSomaticVariantParsing,
                 binding,
                 "Force enable parsing of somatic variants even if the fraction of germline variants is low. Only "
                 "applicable for input VCF that contains flagged somatic variants.")
      ->needs(opt::kVcf)
      ->default_val(false)
      ->group(group_name);
}

// AugmentIntervals Inputs
CLI::Option* AddMappabilityBigwigOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    fs::path& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  // WGS option setup may replace this generic description with library-dependent default details.
  return app
      ->add_option(
          opt::kMappabilityBigwig, binding, "Path to an input BigWig file describing mappability across the genome.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddBlocklistBedOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<fs::path>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kBlocklistBed,
                   binding,
                   "Path to a BED file of 0-based regions to exclude from the analysis before generating intervals. "
                   "For WGS workflows, an omitted value selects a command-specific bundled hg38 blocklist.")
      ->check(cli::BedFileValidator())
      ->group(group_name);
}

// AugmentIntervals Other
CLI::Option* AddIncludeAltContigsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kIncludeAltContigs,
                 binding,
                 "include alt contigs when generating the interval file (default: ignore alt and decoy contigs)")
      ->default_val(kAugmentIntervalsDefaultIncludeAltContigs)
      ->group(group_name);
}

// AugmentIntervals Whole Genome
CLI::Option* AddWholeGenomeOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kWholeGenome,
                 binding,
                 "generate intervals of `--interval-size` across the whole genome. Using this option will "
                 "turn off off-target generation.")
      ->default_val(kAugmentIntervalsDefaultWholeGenome)
      ->group(group_name);
}

CLI::Option* AddWholeGenomeIntervalSizeOption(CLI::App* const app,
                                              const std::string& group_name,
                                              size_t& binding,
                                              const size_t default_size,
                                              const std::string& description_suffix) {
  std::string description =
      "The size (in base pairs) for each target region of the genome for which log ratios are calculated.";
  if (!description_suffix.empty()) {
    description += " " + description_suffix;
  }
  return app->add_option(opt::kIntervalSize, binding, description)
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(default_size)
      ->group(group_name);
}

CLI::Option* AddWholeGenomeIntervalSizeOption(CLI::App* const app,
                                              const std::string& group_name,
                                              size_t& binding,
                                              const size_t default_size) {
  return AddWholeGenomeIntervalSizeOption(app, group_name, binding, default_size, "");
}

CLI::Option* AddMinGcContentOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMinGcContent,
                   binding,
                   "Minimum GC content threshold for filtering targets. "
                   "Targets with GC content below this value will be removed.")
      ->default_val(kAugmentIntervalsDefaultMinGcContent)
      ->check(CLI::Range(0.0, 1.0))
      ->group(group_name);
}

CLI::Option* AddMaxGcContentOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMaxGcContent,
                   binding,
                   "Maximum GC content threshold for filtering targets. "
                   "Targets with GC content above this value will be removed.")
      ->default_val(kAugmentIntervalsDefaultMaxGcContent)
      ->check(CLI::Range(0.0, 1.0))
      ->group(group_name);
}

CLI::Option* AddWholeGenomeMinMappabilityOption(CLI::App* const app,
                                                const std::string& group_name,
                                                f64& binding,
                                                const f64 default_mappability) {
  // non-const required
  return app
      ->add_option(opt::kMinIntervalMappability,
                   binding,
                   "Minimum allowable mappability score for an interval. Intervals with a score less than this "
                   "value are excluded to ensure high-quality read alignment. Mappability scores come from "
                   "`--mappability-bigwig`.")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(default_mappability)
      ->group(group_name);
}

// AugmentIntervals Target Options
CLI::Option* AddBedOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<fs::path>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kBed, binding, "intervals file to augment, in BED format")
      ->check(cli::BedFileValidator())
      ->group(group_name);
}

CLI::Option* AddNoOffTargetsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_flag(opt::kNoOffTargets, binding, "do not generate off-target intervals")
      ->default_val(kAugmentIntervalsDefaultNoOffTargets)
      ->group(group_name);
}

CLI::Option* AddOnTargetIntervalSizeOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kOnTargetIntervalSize,
                   binding,
                   "desired on-target interval size. Large on-target intervals will be split up into sub-interval of "
                   "about this size")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kAugmentIntervalsDefaultOnTargetIntervalSize)
      ->group(group_name);
}

CLI::Option* AddOffTargetMinWidthOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kOffTargetMinWidth, binding, "minimum width for off-targets")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kAugmentIntervalsDefaultOffTargetMinWidth)
      ->group(group_name);
}

CLI::Option* AddOffTargetTrimLengthOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kOffTargetTrimLength,
                   binding,
                   "number of bases to trim off from candidate off-targets to prevent reads overlapping 2+ off-targets")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kAugmentIntervalsDefaultOffTargetTrimLength)
      ->group(group_name);
}

CLI::Option* AddOffTargetIntervalSizeOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kOffTargetIntervalSize,
                   binding,
                   "desired off-target interval size. Large off-target intervals will be split up into sub-intervals "
                   "of about this size, while off-targets smaller than this size will be left as-is (if they exceed "
                   "off-target-min-width)")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kAugmentIntervalsDefaultOffTargetIntervalSize)
      ->group(group_name);
}

CLI::Option* AddOnTargetMinMappabilityOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kOnTargetMinMappability, binding, "minimum mappability threshold for on-target intervals")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kAugmentIntervalsDefaultOnTargetMinMappability)
      ->group(group_name);
}

CLI::Option* AddOffTargetMinMappabilityOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kOffTargetMinMappability, binding, "minimum mappability threshold for off-target intervals")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kAugmentIntervalsDefaultOffTargetMinMappability)
      ->group(group_name);
}

CLI::Option* AddSexChromosomeTargetMinMappabilityOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kSexChromosomeTargetMinMappability,
                   binding,
                   "minimum mappability threshold for sex-chromosome-target intervals. This is will be applied to "
                   "both on- and off-target intervals on sex chromosomes. Thus it overrides the "
                   "--on-target-min-mappability and --off-target-min-mappability values for targets on the sex "
                   "chromosome")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kAugmentIntervalsDefaultSexChromosomeTargetMinMappability)
      ->group(group_name);
}

// CalculateCoverage
CLI::Option* AddCoverageExcludeFlagsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::string& binding) {  // NOSONAR cpp:S995 cpp:M23_131

  return app
      ->add_option(opt::kExcludeFlags,
                   binding,
                   "Exclude alignments if any bits in their FLAG field match the specified integer. Refer to "
                   "https://broadinstitute.github.io/picard/explain-flags.html for flag generation.")
      ->default_val(kCalculateCoverageDefaultExcludeFlags)
      ->group(group_name);
}

CLI::Option* AddCoverageIgnoreDnOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kSkipDelAndRefBases, binding, "Whether to ignore deletion (D) and ref skip (N) in base counting.")
      ->default_val(kCalculateCoverageDefaultIgnoreDN)
      ->group(group_name);
}

// GCCorrect
CLI::Option* AddGcCorrectFirstSpanOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kFirstPassGcSmoothingFraction,
                   binding,
                   "Proportion of data points used for smoothing during the first LOESS iteration of GC bias "
                   "correction.")
      ->check(CLI::TypeValidator<f64>())
      ->default_val(kGCCorrectDefaultFirstSpan)
      ->group(group_name);
}

// TwoSampleLogR / Shared coverage
CLI::Option* AddTumorMinCoverageOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMinTumorCoverage,
                   binding,
                   "Minimum allowable coverage in the tumor sample for an interval to be kept. Intervals below "
                   "this value are excluded.")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kCopyNumberCallerDefaultTumorMinCoverage)
      ->group(group_name);
}

CLI::Option* AddNormalMinCoverageOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMinNormalCoverage,
                   binding,
                   "Minimum allowable coverage in the normal sample for an interval to be kept. Intervals below "
                   "this value are excluded.")
      ->check(CLI::TypeValidator<size_t>())
      ->default_val(kCopyNumberCallerDefaultSomaticNormalMinCoverage)
      ->group(group_name);
}

CLI::Option* AddTumorCoverageFileOption(
    CLI::App* const app,
    const std::string& group_name,
    std::optional<fs::path>&

        // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
        binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kTumorCoverageFile, binding, "tumor coverage file")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddNormalCoverageFileOption(
    CLI::App* const app,
    const std::string& group_name,
    std::optional<fs::path>&

        // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
        binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kNormalCoverageFile, binding, "normal coverage file")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddCoverageFileOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<fs::path>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kCoverageFile, binding, "coverage file")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

// Segmentation
CLI::Option* AddSeedSegmentsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    fs::path& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kSeedSegments,
                   binding,
                   "Path to pre-generated segments from which to start segmentation. The seed segments file must be in "
                   "SEG format. WGS workflows use the bundled hg38 resource at "
                   "/resources/hg38_seed.seg when omitted.")
      ->check(cli::FileReadableValidator())
      ->group(group_name);
}

CLI::Option* AddDisableHierarchicalPruningOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kSkipHierarchicalPruning, binding, "Disable the post-segmentation hierarchical pruning step.")
      ->group(group_name);
}

CLI::Option* AddDisableMergingOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kSkipTTestMerging,
                 binding,
                 "If enabled, the t-test comparison that merges sub-segments back together will be skipped. "
                 "Segments will not have a chance to be re-merged before continuing the recursive segmentation "
                 "process.")
      ->group(group_name);
}

CLI::Option* AddCbsMethodOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    CbsMaxTMethod& binding) {  // NOSONAR cpp:S995 cpp:M23_131
                               // CLI11 on parse; non-const required
  return app
      ->add_option(opt::kCbsMethod,
                   binding,
                   "disable merging at the end of each recursive step. Only applicable if "
                   "--segmentation-mode=germline")
      ->transform(CLI::Transformer(kStringToCbsMethod, CLI::ignore_case))
      ->group(group_name);
}

CLI::Option* AddSegmentationModeOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    SegmentationMode& binding) {  // NOSONAR cpp:S995 cpp:M23_131

  return app->add_option(opt::kSegmentationMode, binding, "choose between germline or somatic")
      ->transform(CLI::Transformer(kStringToSegmentationMode, CLI::ignore_case))
      ->group(group_name);
}

CLI::Option* AddTruncateOutliersOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_flag(opt::kTruncateOutliers, binding, "disable merging at the end of each recursive step")
      ->group(group_name);
}

CLI::Option* AddPruningClusteringParameterOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kHierarchicalPruningHeight, binding, "Clustering parameter for hierarchical pruning.")
      ->default_val(kSegmentationDefaultPruningClusteringParameter)
      ->group(group_name);
}

CLI::Option* AddBafDistanceWeightOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(
          opt::kBafDistanceWeight,
          binding,
          "Weight applied to the mirrored-BAF axis of the distance used for hierarchical pruning, applied only "
          "when both segments in a pair are copy-neutral. Values >1 increase the separation of "
          "allelic-imbalance segments (e.g. copy-neutral LOH) from balanced segments, making them less likely "
          "to be merged away.")
      ->default_val(kSegmentationDefaultBafDistanceWeight)
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

CLI::Option* AddMinTForAutomaticSegmentationOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMinTValueToSkipPermutationTesting,
                   binding,
                   "Minimum t-test value threshold to skip significance testing for segmentation candidates.")
      ->default_val(kSegmentationDefaultMinTForAutomaticSegmentation)
      ->group(group_name);
}

CLI::Option* AddUndoLogrSegmentsOption(CLI::App* const app,
                                       const std::string& group_name,
                                       SegmentUndoMode& binding,
                                       const SegmentUndoMode default_value) {
  return cli::AddEnumOption(
             app,
             opt::kSegmentUndoing,
             binding,
             "Strategy for undoing log-ratio segments (re-merging segments with similar median log-ratio values). "
             "`disabled`: no undoing. `iterative`: perform at least one merge pass at "
             "`--initial-segment-undoing-factor`, then relax the factor each subsequent round (adding "
             "`--segment-undoing-increment` per round), re-merging until `--max-segments-for-undoing` is "
             "reached or no more segments can be undone. The first pass always runs even when the segment count is "
             "already at/below `--max-segments-for-undoing`. `non-iterative`: perform a "
             "single merge pass at `--initial-segment-undoing-factor` (`--max-segments-for-undoing` and "
             "`--segment-undoing-increment` are ignored).",
             default_value)
      ->group(group_name);
}

CLI::Option* AddUndoLogrSegmentsSdFactorOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kInitialSegmentUndoingFactor,
                   binding,
                   "The initial standard deviation multiplier for segment undoing. Adjacent segments with a "
                   "difference less than or equal to the threshold determined by this multiplier are merged. "
                   "Only used when segment undoing is enabled.")
      ->default_val(kSegmentationDefaultUndoLogrSegmentsSdFactor)
      ->check(CLI::TypeValidator<f64>())
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

CLI::Option* AddIncrementUndoLogrSegmentsSdFactorOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(
          opt::kSegmentUndoingIncrement,
          binding,
          "Step by which `--initial-segment-undoing-factor` is relaxed during `iterative` undoing. Each round after "
          "the first adds this value to the factor, gradually relaxing the merging criteria until the stop condition "
          "is met. Only used by the `iterative` strategy.")
      ->default_val(kSegmentationDefaultIncrementUndoLogrSegmentsSdFactor)
      ->check(CLI::TypeValidator<f64>())
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

CLI::Option* AddMaxNumSegmentsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMaxSegmentsForUndoing,
                   binding,
                   "Target segment count for the segment undoing procedure. The process stops once the total "
                   "number of segments is less than or equal to this value. "
                   "Only used when segment undoing is enabled.")
      ->default_val(kSegmentationDefaultMaxNumSegments)
      ->check(CLI::TypeValidator<size_t>())
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

CLI::Option* AddUndoDhSegmentsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kEnableDhSegmentUndoing,
                 binding,
                 "Undo segments based on DH values (re-merge segments with similar median DH values).")
      ->default_val(kSegmentationDefaultUndoDhSegments)
      ->group(group_name);
}

CLI::Option* AddUndoDhSegmentsSdFactorOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kInitialDhSegmentUndoingFactor,
                   binding,
                   "The initial standard deviation multiplier for DH segment undoing. Adjacent segments with a "
                   "difference less than or equal to the threshold determined by this multiplier are merged. "
                   "Requires `--enable-dh-segment-undoing`.")
      ->default_val(kSegmentationDefaultUndoDhSegmentsSdFactor)
      ->group(group_name);
}

CLI::Option* AddMinObservationsPerSegmentOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(
          opt::kMinObservationsPerSegment, binding, "Minimum number of observations a segment is allowed to have")
      ->default_val(kSegmentationDefaultMinObsPerSegment)
      ->group(group_name);
}

// PurityPloidySearch
CLI::Option* AddSegMaxLogrOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMaxSegmentLogRatioForPurityPloidySearch,
                   binding,
                   "Maximum log-ratio value allowed for a segment to be included in the purity/ploidy search. "
                   "Segments with a value greater than this are excluded.")
      ->default_val(kSegMaxLogR)
      ->group(group_name);
}

CLI::Option* AddSegMinNumLogrsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMinLogRatiosPerSegmentForPurityPloidySearch,
                   binding,
                   "Minimum number of log-ratio observations required to include a segment in the purity/ploidy "
                   "search. Segments with fewer observations than this value are excluded.")
      ->default_val(kSegMinNumLogRs)
      ->group(group_name);
}

CLI::Option* AddSegMinNumSnpsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMinSnpsPerSegmentForPurityPloidySearch,
                   binding,
                   "Minimum number of SNPs required to consider a segment in the purity/ploidy search. Segments "
                   "with fewer SNPs than this value are excluded.")
      ->default_val(kSegMinNumSnps)
      ->group(group_name);
}

CLI::Option* AddBafLikelihoodWeightForPurityPloidySearchRerunOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding,
    const f64 default_value) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kBafLikelihoodWeightForPurityPloidySearchRerun,
                   binding,
                   "Weight applied to the BAF term of the joint likelihood (joint = log-ratio + weight * BAF) "
                   "during the purity-ploidy-search re-run only. The initial search is unaffected. "
                   "A value greater than 1.0 trusts the BAF over the log ratio")
      ->check(CLI::TypeValidator<f64>())
      ->check(CLI::NonNegativeNumber)
      ->default_val(default_value)
      ->group(group_name);
}

CLI::Option* AddMinLohFractionForWgdOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMinLohFractionForWgd,
                   binding,
                   "Minimum fraction of the assessable autosomal genome (by length) that must be "
                   "loss-of-heterozygosity for the low-purity whole-genome-doubling "
                   "(WGD) "
                   "fallback to reject the diploid-bias solution. Only used when tumor purity is too low for the "
                   "mBAF-peak WGD test")
      ->check(CLI::TypeValidator<f64>())
      ->check(CLI::Range(0.0, 1.0))
      ->default_val(kMinLohFractionForWgd)
      ->group(group_name);
}

CLI::Option* AddTumorOnlyOption(CLI::App* const app,
                                const std::string& group_name,
                                bool& binding,
                                const std::string& description) {
  return app->add_flag(opt::kTumorOnly, binding, description)->default_val(false)->group(group_name);
}

// Sample Metadata
CLI::Option* AddSexOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<sex_predict::Sex>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  const auto desc = fmt::format("value in {}", cli::FormatEnumNames<sex_predict::Sex>("{", "}"));
  return app
      ->add_option_function<std::string>(
          opt::kSex,
          [&binding](const std::string& input) {
            binding = cli::ParseEnumNameOrThrow<sex_predict::Sex>(std::string{opt::kSex}, input);
          },
          "Sex of sample.")
      ->transform([](const std::string& value) -> std::string { return value; }, desc)
      ->group(group_name);
}

CLI::Option* AddSexSelectionOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    sex_predict::SexSelection& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  const auto desc = fmt::format("value in {}, or deprecated alias unknown",
                                cli::FormatEnumNames<sex_predict::SexSelection>("{", "}"));
  return app
      ->add_option_function<std::string>(
          opt::kSex,
          [&binding](const std::string& input) {
            const auto maybe_sex = enum_util::ParseEnumName<sex_predict::Sex>(input);
            if (maybe_sex.has_value() && maybe_sex.value() == sex_predict::Sex::kUnknown) {
              Logging::Warn("'--sex unknown' is deprecated; use '--sex skip-prediction'");
              binding = sex_predict::SexSelection::kSkipPrediction;
              return;
            }
            binding = cli::ParseEnumNameOrThrow<sex_predict::SexSelection>(std::string{opt::kSex}, input);
          },
          "How to obtain the sample sex.")
      ->transform([](const std::string& value) { return value; }, desc)
      ->default_val(cli::FormatEnumName(sex_predict::SexSelection::kAutoDetect))
      ->group(group_name);
}

CLI::Option* AddLibraryTypeOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<io::LibraryType>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  const auto desc = fmt::format("value in {}", cli::FormatEnumNames<io::LibraryType>("{", "}"));
  return app
      ->add_option_function<std::string>(
          opt::kLibraryType,
          [&binding](const std::string& input) {
            binding = cli::ParseEnumNameOrThrow<io::LibraryType>(std::string{opt::kLibraryType}, input);
          },
          "Sequencing library type of the input BAM. Used to select the default interval size and nonIntegerTcn "
          "thresholds for germline WGS. When unset, it is auto-detected from the BAM.")
      ->transform([](const std::string& value) { return value; }, desc)
      ->group(group_name);
}

CLI::Option* AddSampleIdOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::string& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kSampleId, binding, "Sample ID.")
      ->check(CLI::TypeValidator<std::string>())
      ->group(group_name);
}

CLI::Option* AddSampleNameOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<std::string>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(
          opt::kSampleName,
          binding,
          "Name of sample in the input VCF header and used as a sample ID in the output. Does not require `--vcf`.")
      ->check(CLI::TypeValidator<std::string>())
      ->group(group_name);
}

CLI::Option* AddNormalSampleNameOption(
    CLI::App* const app,
    const std::string& group_name,
    std::optional<std::string>&

        // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
        binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kNormalSampleName, binding, "Normal sample name.")
      ->check(CLI::TypeValidator<std::string>())
      ->group(group_name);
}

CLI::Option* AddTumorSampleNameOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<std::string>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kTumorSampleName, binding, "Tumor sample name.")
      ->check(CLI::TypeValidator<std::string>())
      ->group(group_name);
}

CLI::Option* AddTumorPurityOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<f64>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  constexpr f64 kMinTumorPurity = 0.15;
  return app->add_option(opt::kTumorPurity, binding, "tumor purity")
      ->check(CLI::TypeValidator<f64>())
      ->check(CLI::Range(static_cast<f64>(kMinTumorPurity), static_cast<f64>(1)))
      ->group(group_name);
}

CLI::Option* AddTumorPloidyOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::optional<f64>& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kTumorPloidy, binding, "tumor ploidy")
      ->check(CLI::TypeValidator<f64>())
      ->check(CLI::PositiveNumber)
      ->group(group_name);
}

// Likelihood
CLI::Option* AddCnvLengthFlagMinSizeOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    size_t& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kMinCnvLengthToFlag, binding, "Flag calls with a size lower than this threshold (2000).")
      ->default_val(kLikelihoodCnvLengthFlagMinSize)
      ->group(group_name);
}

// Save-flag helpers (new hidden advanced options)
CLI::Option* AddSaveTaskflowGraphOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kSaveTaskflowGraph, binding, "Save intermediate taskflow graph output helpful for debugging.")
      ->group(group_name);
}

CLI::Option* AddSaveLogRatioSegmentsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kSaveLogRatioSegments,
                 binding,
                 "Output file for pre-likelihood LogR segments. These segments have gone through the undoing and "
                 "pruning processes already, but have not merged based on equal copy number or by mapping "
                 "qualities.")
      ->group(group_name);
}

CLI::Option* AddSaveBafSegmentsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kSaveBafSegments,
                 binding,
                 "Output file for segments after BAF segmentation but before hierarchical pruning.")
      ->group(group_name);
}

CLI::Option* AddWriteSegmentLikelihoodOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kWriteSegmentLikelihood,
                 binding,
                 "Emit per-segment, per-purity/ploidy-grid-cell likelihood diagnostics, one TSV per search pass "
                 "(segment_likelihood.diploid_bias.tsv and, when the WGD QC re-run triggers, "
                 "segment_likelihood.uniform.tsv).")
      ->group(group_name);
}

CLI::Option* AddMapqCutoffToFlagCallsOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    s32& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kMinMapqForCalls,
                   binding,
                   "Minimum average mapping quality (MAPQ) required for a call in the SEG output. Segments with "
                   "an average MAPQ less than this value are merged. If merging is not possible, the call will be "
                   "flagged in the SEG output. Use a negative value to disable this filter.")
      ->default_val(kLikelihoodDefaultGermlineMAPQCutoffForCalls)
      ->group(group_name);
}

// MergeSegments
CLI::Option* AddMinMapqOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    f64& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(
          opt::kMinMapq, binding, "minimum avg mean MAPQ for keeping segments. values less than this will get merged")
      ->default_val(kMergeDefaultMinMapq)
      ->group(group_name);
}

CLI::Option* AddUseMapqsFromFileOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kUseMapqsFromFile,
                 binding,
                 "use the mean mapqs per target from the `--mapqs` file instead of the avg mean MAPQs per segment "
                 "contained the segments file. Can use this if the segment file doesn't contain a meanAvgMapq field")
      ->needs(opt::kMapqs)
      ->default_val(kMergeSegmentsDefaultUseMapqsObservations)
      ->group(group_name);
}

CLI::Option* AddRecalculatePerSegmentDataOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_flag(opt::kRecalculatePerSegmentData,
                 binding,
                 "recalculate the following fields from the SEG file: NumLogRatio, MeanLogRatio, NumSnp, MeanDh, "
                 "AvgMeanMapq etc. Otherwise, calculates the weighted mean of these fields for each merged segment.")
      ->needs(opt::kLogRatios)
      ->needs(opt::kMapqs)
      ->needs(opt::kSex)
      ->default_val(kMergeSegmentsDefaultRecalculatePerSegmentData)
      ->group(group_name);
}

CLI::Option* AddOutputDirOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    fs::path& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_option(opt::kOutputDir, binding, "Path to output directory.")->default_val(".")->group(group_name);
}

CLI::Option* AddOutputPrefixOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_option stores a non-const reference for later population during parse; opts must remain mutable.
    std::string& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app
      ->add_option(opt::kOutputPrefix,
                   binding,
                   "Prefix prepended to all output filenames (e.g. a sample ID), joined with an underscore. "
                   "If unset, default filenames are used.")
      ->check(cli::FilenameComponentValidator())
      ->group(group_name);
}

CLI::Option* AddOverwriteOption(
    CLI::App* const app,
    const std::string& group_name,
    // CLI11 add_flag stores a non-const reference for later population during parse; opts must remain mutable.
    bool& binding) {  // NOSONAR cpp:S995 cpp:M23_131
  return app->add_flag(opt::kOverwrite, binding, "Overwrite output files if they already exist.")->group(group_name);
}

}  // namespace xoos::cnc
