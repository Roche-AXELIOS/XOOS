#include "demux.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <optional>
#include <string>

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/file-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/nonempty-file-validator.h>
#include <xoos/enum/enum-util.h>
#include <xoos/util/file-functions.h>
#include <xoos/vfs/vfs.h>

#include "adapter-design/adapter-design-constants.h"
#include "adapter-dsl/authored/authored-json.h"
#include "io/sequence-reader.h"
#include "task/batch.h"
#include "task/flow-context.h"

namespace xoos::demux {

using cli::AddEnumOption;
using cli::AddOptionalEnumOption;
using cli::AddOutputFileOption;
using cli::AddThreadCountOption;

namespace {

// Store deprecated values separately so the replacement option takes precedence when both names are specified.
struct DeprecatedDemuxCliOptions {
  xoos::vec<fs::path> reads;
  xoos::vec<fs::path> short_reads;
  std::optional<fs::path> reads_file_list;
  std::optional<fs::path> output_dir;
  std::optional<fs::path> short_output_dir;
  std::optional<size_t> writers_per_sample;
  std::optional<size_t> workers_per_reads_file;
  std::optional<size_t> min_trimmed_read_length;
  bool write_failed_reads{false};
};

using DeprecatedDemuxCliOptionsPtr = std::shared_ptr<DeprecatedDemuxCliOptions>;

// Register a hidden removed option that reports its replacement when used.
void AddRemovedOption(CLI::App* const app, const std::string& old_name, const std::string& replacement) {
  app->add_option_function<std::string>(
         old_name,
         [old_name, replacement](const std::string&) {
           throw CLI::ValidationError(
               fmt::format("The option '{}' is no longer supported; use '{}' instead.", old_name, replacement));
         },
         "")
      ->group("");
}

// Warn when a deprecated option is used and when its replacement takes precedence.
void WarnDeprecatedOption(const CLI::App* const app,
                          const std::string& old_name,
                          const std::string& replacement,
                          const bool replacement_specified,
                          const std::string& conflicting_option) {
  if (app->count(old_name) == 0) {
    return;
  }
  if (replacement_specified || app->count(replacement) > 0) {
    Logging::Warn("The deprecated option '{}' is ignored because '{}' was also specified; use '{}' only.",
                  old_name,
                  conflicting_option.empty() ? replacement : conflicting_option,
                  replacement);
    return;
  }
  Logging::Warn("The option '{}' is deprecated and will be removed in a future release; use '{}' instead.",
                old_name,
                replacement);
}

// Use this overload when only the deprecated option is relevant and no precedence conflict needs reporting.
void WarnDeprecatedOption(const CLI::App* const app, const std::string& old_name, const std::string& replacement) {
  WarnDeprecatedOption(app, old_name, replacement, false, {});
}

// Apply a deprecated reads option after expanding paths, unless a preferred reads option was specified.
void ApplyDeprecatedReads(const CLI::App* const app,
                          const std::string& old_name,
                          const std::string& replacement,
                          const vec<fs::path>& deprecated_reads,
                          const bool replacement_specified,
                          const std::string& conflicting_option,
                          DemuxAndTrimParam& param) {
  WarnDeprecatedOption(app, old_name, replacement, replacement_specified, conflicting_option);
  if (app->count(old_name) == 0 || replacement_specified) {
    return;
  }
  for (const auto& input : deprecated_reads) {
    if (!fs::exists(input)) {
      throw CLI::ValidationError(fmt::format("Input does not exist: '{}'", input.string()));
    }
  }
  const auto inputs = ExpandInputFileList(deprecated_reads);
  std::ranges::copy(inputs, std::back_inserter(param.inputs));
}

// Apply a deprecated option's value only when its replacement was not specified.
template <typename T>
void ApplyDeprecatedValue(const CLI::App* const app,
                          const std::string& old_name,
                          const std::string& replacement,
                          const std::optional<T>& deprecated_value,
                          T& value) {
  WarnDeprecatedOption(app, old_name, replacement);
  if (app->count(old_name) > 0 && app->count(replacement) == 0) {
    value = deprecated_value.value();
  }
}

}  // namespace

// Return all sequence files found by recursively traversing a path.
static vec<fs::path> RecursiveFindSequenceFiles(const fs::path& input) {
  vec<fs::path> stack{input};

  vec<fs::path> files;
  while (!stack.empty()) {
    const auto current = stack.back();
    stack.pop_back();

    if (fs::is_directory(current)) {
      for (const auto& file : fs::directory_iterator(current)) {
        stack.push_back(file.path());
      }
    } else if (IsSequenceFileFormat(current)) {
      files.push_back(current);
    }
  }
  return files;
}

static vec<fs::path> ReadInputFileList(const fs::path& input_file_list) {
  std::ifstream ifs{input_file_list};
  vec<fs::path> input_files;
  std::string line;
  while (std::getline(ifs, line)) {
    if (!fs::exists(line)) {
      throw CLI::ValidationError(fmt::format("Input does not exist: '{}'", line));
    }
    input_files.emplace_back(line);
  }
  return input_files;
}

vec<fs::path> ExpandInputFileList(const vec<fs::path>& input_files) {
  vec<fs::path> expanded_input_files;
  for (const auto& input_file : input_files) {
    // This public helper can be called without the CLI::ExistingPath validator, so ignore missing paths here.
    if (!fs::exists(input_file)) {
      continue;
    }
    // Do not apply directory filename filtering to explicit files. Keeping the path lets OpenSequenceFile report an
    // unsupported format instead of silently dropping the input; directory discovery still ignores unrelated files.
    if (!fs::is_directory(input_file)) {
      expanded_input_files.push_back(input_file);
      continue;
    }
    const auto files{RecursiveFindSequenceFiles(input_file)};
    expanded_input_files.insert(expanded_input_files.end(), files.cbegin(), files.cend());
  }
  return expanded_input_files;
}

// Validate --adapter-design-name and enable methylation mode when SBX-DM is selected.
static std::string ValidateAdapterDesignName(DemuxAndTrimParam& param, const std::string_view input) {
  if (input.empty()) {
    return "adapter-design-name must not be empty";
  }
  if (input == "SBX-DM") {
    param.methylation = true;
  }
  return {};
}

// Read the adapter names from the selected bundle so a missing required value has an actionable error message.
static std::string FormatAdapterDesignNames(const fs::path& adapter_design_bundle) {
  const auto bundle_vfs = vfs::Open(adapter_design_bundle);
  const dsl::AdapterBundle bundle = dsl::FromAuthoredJsonBundle(bundle_vfs, fs::path{kAdapterConfigFilename});
  vec<std::string> names;
  names.reserve(bundle.adapters.size());
  for (const auto& [name, _] : bundle.adapters) {
    names.push_back(name);
  }
  std::ranges::sort(names, std::greater<>{});
  return fmt::format("{}", fmt::join(names, ", "));
}

// Configure strand detection from a precomputed index using its average FPR and the configured critical value.
static void ConfigureStrandDetector(DemuxAndTrimParam& param, const fs::path& filename) {
  param.strand_detector.emplace(filename);
  const auto [fw_fpr, rv_fpr] = param.strand_detector->FalsePositiveRates();
  Logging::Info("Strand Detection Mode. Filter FPRs fw: {:.5f} rv: {:.5f}", fw_fpr, rv_fpr);
  const f64 average_error_rate = (fw_fpr + rv_fpr) / 2.0;
  const f64 rescaled_critical_value = std::pow(10, -param.strand_critical / 10.0);
  param.strand_classifier.emplace(*param.strand_detector, average_error_rate, rescaled_critical_value);
}

// Define CLI options specific to simplex adapter processing.
static void DefineSimplexOptions(CLI::App* const app, const std::shared_ptr<DemuxAndTrimParam>& param) {
  AddEnumOption(app,
                "--read-length-mode",
                param->read_length_mode,
                fmt::format("Simplex read output mode. "
                            "{}: emit full and partial reads together. "
                            "{}: emit only full-length reads. "
                            "{}: emit full and partial reads into separate subdirectories (<sample>/{}/, <sample>/{}).",
                            enum_util::FormatEnumName(ReadLengthMode::kAll),
                            enum_util::FormatEnumName(ReadLengthMode::kFullOnly),
                            enum_util::FormatEnumName(ReadLengthMode::kAllSplit),
                            kFullReadSubdir,
                            kPartialReadSubdir),
                param->read_length_mode)
      ->group("Simplex Options");
  AddEnumOption(app,
                "--discordant-sid-mode",
                param->discordant_sid_mode,
                fmt::format("Determines how to handle reads where 5' and 3' SIDs disagree. "
                            "{}: discard only when both sides are equally likely (tied). "
                            "{}: discard all discordant reads (both sides significantly detected). "
                            "{}: keep all discordant reads, assign to winning side. Ties are set to 5'.",
                            enum_util::FormatEnumName(DiscordantSidMode::kDiscardTied),
                            enum_util::FormatEnumName(DiscordantSidMode::kDiscardAll),
                            enum_util::FormatEnumName(DiscordantSidMode::kKeep)),
                param->discordant_sid_mode)
      ->group("Simplex Options");
  app->add_option_function<std::string>(
         "--dual-sid-mode",
         [&param](const std::string& value) {
           const auto parsed = ParseDualSidMode(value);
           if (!parsed) {
             throw CLI::ValidationError(
                 fmt::format("dual-sid-mode: Check {} value in {{none, need-5p, need-3p, need-both}} FAILED", value));
           }
           param->dual_sid_mode = *parsed;
         },
         fmt::format("Controls single-ended assignment for dual-SID simplex sample sheets. "
                     "{}: derive 5p-only and 3p-only eligibility from per-end SID collisions. "
                     "{}: require the 5p SID side for assignment. "
                     "{}: require the 3p SID side for assignment. "
                     "{}: require both SID sides for assignment.",
                     ToString(DualSidMode::kNone),
                     ToString(DualSidMode::kNeed5p),
                     ToString(DualSidMode::kNeed3p),
                     ToString(DualSidMode::kNeedBoth)))
      ->default_val(std::string(ToString(param->dual_sid_mode)))
      ->group("Simplex Options");
  app->add_flag("--override-simplex-quality",
                param->override_simplex_qual_score,
                "Override base-call quality scores for simplex adapter reads (YS, YSU, YS-NEW, SIMPLEX-10X) "
                "with configured simplex value (see --quality-score-overrides; defaults to '7'). "
                "By default, the original base-call quality scores are preserved.")
      ->default_val(param->override_simplex_qual_score)
      ->group("Simplex Options");
  app->add_option("--min-score", param->min_score, "Minimum log odds based score needed for a valid adapter match.")
      ->check(CLI::Range(s32{0}, std::numeric_limits<s32>::max()))
      ->default_val(param->min_score)
      ->group("Simplex Options");
}

// Define CLI options specific to duplex adapter processing and strand detection.
static void DefineDuplexOptions(CLI::App* const app, const std::shared_ptr<DemuxAndTrimParam>& param) {
  app->add_option("--max-pairwise-error-rate",
                  param->max_pairwise_error_rate,
                  "Maximum pairwise error rate for retaining a duplex consensus region.")
      ->check(CLI::Range(0.0F, 1.0F))
      ->default_val(param->max_pairwise_error_rate)
      ->group("Duplex Options");
  AddRemovedOption(app, "--max-error-rate-percent", "--max-pairwise-error-rate");
  app->add_option("--min-concordant-duplex-bases-for-early-stop",
                  param->min_concord_dp_bases,
                  "Minimum concordant duplex bases per sample for early stopping.")
      ->group("Duplex Options");
  AddRemovedOption(app, "--stop-after-min-concordant-duplex-bases", "--min-concordant-duplex-bases-for-early-stop");
  app->add_option("--strand-critical-phred",
                  param->strand_critical,
                  "In duplex strand detection mode, this is the Phred scaled FPR critical value "
                  "-10*log10(critical_value). Low values will classify faster at cost of accuracy.")
      ->default_val(param->strand_critical)
      ->check(CLI::PositiveNumber)
      ->group("Duplex Options");
  app->add_option_function<fs::path>(
         "--strand-index",
         [&param](const fs::path& filename) { ConfigureStrandDetector(*param, filename); },
         "Enable genome strand detection with the specified file created by demux_strand_index.")
      ->check(cli::NonEmptyFileValidator())
      ->group("Duplex Options");
  AddRemovedOption(app, "--strand-detect", "--strand-index");
  AddEnumOption(app,
                "--ambiguous-topology-quality-mode",
                param->ambiguous_topology_quality_mode,
                fmt::format("Mode to set base-call quality score (default: '{}') to duplex consensus bases at "
                            "topologically ambiguous indel positions. "
                            "{}: disabled; ambiguous bases keep their normal duplex quality score. "
                            "{}: smear the ambiguous score over concordant bases that an indel could validly occupy "
                            "(consensus-preserving; populates the ambiguous_duplex_bases metric). "
                            "{}: {} plus rewrite discordant bases to the ambiguous score where the smear marks them "
                            "(overridden bases count as ambiguous; the YC tag is unchanged).",
                            BaseQualityScores::kAmbiguousBaseQual,
                            enum_util::FormatEnumName(AmbiguousTopologyQualityMode::kNone),
                            enum_util::FormatEnumName(AmbiguousTopologyQualityMode::kSmear),
                            enum_util::FormatEnumName(AmbiguousTopologyQualityMode::kSmearOverride),
                            enum_util::FormatEnumName(AmbiguousTopologyQualityMode::kSmear)),
                param->ambiguous_topology_quality_mode)
      ->group("Duplex Options");
  app->add_option_function<std::string>(
      "--quality-score-overrides",
      [&param](const std::string& value) {
        if (value.size() != BaseQualityScores::kQualityBaseValues) {
          throw CLI::ValidationError(fmt::format(
              "--quality-score-overrides must be exactly 4 characters (got {}): '{}'", value.size(), value));
        }
        for (const char c : value) {
          if (c < '!' || c > '~') {
            throw CLI::ValidationError(fmt::format(
                "--quality-score-overrides contains a non phred+33 character (must be '!'..'~'): '{}'", value));
          }
        }
        param->base_quality_scores = BaseQualityScores{.simplex = value[BaseQualityScores::kSimplexIndex],
                                                       .concordant = value[BaseQualityScores::kConcordantIndex],
                                                       .discordant = value[BaseQualityScores::kDiscordantIndex],
                                                       .ambiguous = value[BaseQualityScores::kAmbiguousIndex]};
      },
      fmt::format("Output base-call quality scores (phred+33) as a 4-character string in the order "
                  "simplex, concordant duplex, discordant duplex, ambiguous duplex. "
                  "Each character must be in the printable phred+33 range ('!'..'~', Q0..Q93). "
                  "The simplex score applies to both simplex reads and the simplex tail of duplex reads. "
                  "[default '{}{}{}{}']",
                  BaseQualityScores::kSimplexBaseQual,
                  BaseQualityScores::kConcordBaseQual,
                  BaseQualityScores::kDiscordBaseQual,
                  BaseQualityScores::kAmbiguousBaseQual));
}

// Define input and output CLI options shared by simplex and duplex adapter processing.
static void DefineInputOutputOptions(CLI::App* const app,
                                     const std::shared_ptr<DemuxAndTrimParam>& param,
                                     const DeprecatedDemuxCliOptionsPtr& deprecated) {
  app->add_option_function<vec<fs::path>>(
         "--reads",
         [&param](const vec<fs::path>& inputs) {
           const auto expanded{ExpandInputFileList(inputs)};
           std::ranges::copy(expanded, std::back_inserter(param->inputs));
         },
         "Untrimmed FASTQ, FASTQ.gz, or RDB inputs. Directories are expanded recursively. "
         "May be repeated or space separated.")
      ->check(CLI::ExistingPath);
  app->add_option_function<vec<fs::path>>(
         "--input", [&deprecated](const vec<fs::path>& inputs) { deprecated->reads = inputs; }, "")
      ->group("");
  app->add_option_function<vec<fs::path>>(
         "-i", [&deprecated](const vec<fs::path>& inputs) { deprecated->short_reads = inputs; }, "")
      ->group("");
  app->add_option_function<fs::path>(
         "--reads-file-list",
         [&param](const fs::path& input_file_list) {
           const auto inputs = ExpandInputFileList(ReadInputFileList(input_file_list));
           std::ranges::copy(inputs, std::back_inserter(param->inputs));
         },
         "File listing read paths instead of setting them via --reads.")
      ->check(cli::NonEmptyFileValidator());
  app->add_option("--input-files-list", deprecated->reads_file_list, "")->group("");
  app->add_option("--output-dir", param->out_dir, "Output directory")->default_val(param->out_dir);
  app->add_option("--out-dir", deprecated->output_dir, "")->group("");
  app->add_option("-o", deprecated->short_output_dir, "")->group("");
  app->add_flag("--overwrite", param->overwrite, "Overwrites existing non-empty output directory and its contents.")
      ->default_val(param->overwrite);
}

// Define processing configuration CLI options shared by simplex and duplex adapter processing.
static void DefineProcessingOptions(CLI::App* const app,
                                    const std::shared_ptr<DemuxAndTrimParam>& param,
                                    const DeprecatedDemuxCliOptionsPtr& deprecated) {
  app->add_option("--batch-size", param->batch_size, "Number of reads to process in each batch.")
      ->check(CLI::Range(size_t{1}, std::numeric_limits<size_t>::max()))
      ->default_val(param->batch_size);
  AddRemovedOption(app, "-b", "--batch-size");
  AddThreadCountOption(app, "--threads", param->threads);
  app->add_option("--adapter-design-bundle",
                  param->adapter_design_bundle,
                  "ZIP file or directory containing the adapter architecture and search strategy.")
      ->default_val(param->adapter_design_bundle)
      ->check(CLI::ExistingPath);
  AddRemovedOption(app, "-p", "--adapter-design-bundle");
  app->add_option("--adapter-design-name",
                  param->adapter_design_name,
                  "The adapter design name to be loaded from the bundle's DSL config. Required.")
      ->check(CLI::Validator([&param](const std::string& input) { return ValidateAdapterDesignName(*param, input); },
                             "valid-adapter-design-name"));
  AddRemovedOption(app, "-n", "--adapter-design-name");
  app->add_option("--sample-sheet", param->sample_sheet, "Limit adapter SIDs to those in this sample sheet.")
      ->check(cli::NonEmptyFileValidator());
  app->add_option(
         "--compression-level",
         param->compression_level,
         "Output FASTQ compression level: gzip 1-9 or zstd 1-19. Higher levels trade speed for compression; ignored "
         "when compression is disabled. Uses the algorithm default when omitted.")
      ->check(CLI::Range(kMinCompressionLevel, kMaxZstdCompressionLevel));
  AddThreadCountOption(app, "--writers-per-sample", param->writing_threads_per_sample)
      ->description("Number of writing threads used per sample.")
      ->check(CLI::Range(size_t{1}, BatchTask::kMaxNumberOfSinkWorkers))
      ->default_val(param->writing_threads_per_sample);
  app->add_option("--writing-threads-per-sample", deprecated->writers_per_sample, "")
      ->check(CLI::Range(size_t{1}, BatchTask::kMaxNumberOfSinkWorkers))
      ->group("");
  app->add_option("--workers-per-reads-file",
                  param->worker_threads_per_input,
                  "Controls how many input files are processed concurrently: "
                  "concurrent inputs = threads / workers-per-reads-file (minimum 1). "
                  "Higher values reduce concurrency, which can lower memory pressure for large files.")
      ->check(CLI::Range(size_t{1}, std::numeric_limits<size_t>::max()))
      ->default_val(param->worker_threads_per_input);
  app->add_option("--worker-threads-per-input", deprecated->workers_per_reads_file, "")
      ->check(CLI::Range(size_t{1}, std::numeric_limits<size_t>::max()))
      ->group("");
  AddEnumOption(
      app,
      "--compression-type",
      param->compression_type,
      "Compression type for output FASTQs. If set to none, output will ignore the compression level parameter.",
      param->compression_type);
}

// Define the remaining CLI options shared by simplex and duplex adapter processing.
static void DefineGeneralOptions(CLI::App* const app,
                                 const std::shared_ptr<DemuxAndTrimParam>& param,
                                 const DeprecatedDemuxCliOptionsPtr& deprecated) {
  DefineInputOutputOptions(app, param, deprecated);
  DefineProcessingOptions(app, param, deprecated);
  app->add_option("--max-read-length-bin",
                  param->length_distribution_report_max,
                  "The maximum size the largest bucket in the length histogram distribution report can be.")
      ->check(CLI::Range(u32{1}, std::numeric_limits<u32>::max()))
      ->default_val(param->length_distribution_report_max);
  AddRemovedOption(app, "--length-distribution-report-max", "--max-read-length-bin");
  app->add_option("--min-read-length", param->min_read_len, "Filter raw reads shorter than this length.")
      ->check(CLI::Range(size_t{1}, std::numeric_limits<size_t>::max()))
      ->default_val(param->min_read_len);
  AddRemovedOption(app, "--min-read-len", "--min-read-length");
  app->add_option(
         "--min-trimmed-read-length", param->min_trimmed_read_len, "Filter trimmed reads shorter than this length.")
      ->check(CLI::Range(size_t{1}, std::numeric_limits<size_t>::max()))
      ->default_val(param->min_trimmed_read_len);
  app->add_option("--min-trimmed-read-len", deprecated->min_trimmed_read_length, "")
      ->check(CLI::Range(size_t{1}, std::numeric_limits<size_t>::max()))
      ->group("");
  app->add_flag("--write-failed-reads",
                param->output_failed_reads,
                fmt::format("Output unaltered failed reads into the '{}' subdirectory under the output directory. "
                            "Files are named with prefix '{}'.",
                            kReservedRawFailedName,
                            kReservedRawFailedName))
      ->default_val(param->output_failed_reads);
  app->add_flag("--output-failed-reads", deprecated->write_failed_reads, "")->group("");
  // Hidden option for emitting debug tags in output FASTQ comments for trimming evaluation.
  // Accepts comma-separated values, e.g. --trimming-debug raw-seq
  AddOptionalEnumOption(app, "--trimming-debug", param->trimming_debug, "Trimming debug mode")->group("");
}

static void ApplyDeprecatedOptions(const CLI::App* const app,
                                   const std::shared_ptr<DemuxAndTrimParam>& param,
                                   const DeprecatedDemuxCliOptionsPtr& deprecated) {
  const bool reads_specified = app->count("--reads") > 0;
  const bool short_reads_specified = app->count("-i") > 0;
  ApplyDeprecatedReads(app, "-i", "--reads", deprecated->short_reads, reads_specified, {}, *param);
  ApplyDeprecatedReads(app,
                       "--input",
                       "--reads",
                       deprecated->reads,
                       reads_specified || short_reads_specified,
                       short_reads_specified ? "-i" : "",
                       *param);

  WarnDeprecatedOption(app, "--input-files-list", "--reads-file-list");
  if (app->count("--input-files-list") > 0 && app->count("--reads-file-list") == 0) {
    const auto& input_file_list = deprecated->reads_file_list.value();
    if (!fs::exists(input_file_list)) {
      throw CLI::ValidationError(fmt::format("Input does not exist: '{}'", input_file_list.string()));
    }
    const auto inputs = ExpandInputFileList(ReadInputFileList(input_file_list));
    std::ranges::copy(inputs, std::back_inserter(param->inputs));
  }

  ApplyDeprecatedValue(app, "--out-dir", "--output-dir", deprecated->output_dir, param->out_dir);
  ApplyDeprecatedValue(app, "-o", "--output-dir", deprecated->short_output_dir, param->out_dir);
  ApplyDeprecatedValue(app,
                       "--writing-threads-per-sample",
                       "--writers-per-sample",
                       deprecated->writers_per_sample,
                       param->writing_threads_per_sample);
  ApplyDeprecatedValue(app,
                       "--worker-threads-per-input",
                       "--workers-per-reads-file",
                       deprecated->workers_per_reads_file,
                       param->worker_threads_per_input);
  ApplyDeprecatedValue(app,
                       "--min-trimmed-read-len",
                       "--min-trimmed-read-length",
                       deprecated->min_trimmed_read_length,
                       param->min_trimmed_read_len);

  WarnDeprecatedOption(app, "--output-failed-reads", "--write-failed-reads");
  if (app->count("--output-failed-reads") > 0 && app->count("--write-failed-reads") == 0) {
    param->output_failed_reads = deprecated->write_failed_reads;
  }
}

static void ValidateOptions(const std::shared_ptr<DemuxAndTrimParam>& param) {
  if (param->adapter_design_name.empty()) {
    throw CLI::ValidationError(fmt::format("--adapter-design-name is required. Valid values in '{}': {}",
                                           kAdapterConfigFilename,
                                           FormatAdapterDesignNames(param->adapter_design_bundle)));
  }
  // Make sure that we have at least one input file.
  if (param->inputs.empty()) {
    throw CLI::ValidationError("Must provide at least one input");
  }
  if (param->compression_type == SinkCompressionType::kGzip && param->compression_level.has_value() &&
      *param->compression_level > kMaxGzipCompressionLevel) {
    throw CLI::ValidationError(fmt::format("Invalid gzip compression level: {}, must be in the range [{}, {}]",
                                           *param->compression_level,
                                           kMinCompressionLevel,
                                           kMaxGzipCompressionLevel));
  }
  // Sort inputs for deterministic output filename assignment.
  std::ranges::sort(param->inputs);
  // Guard against hardware_concurrency() returning 0 when --threads 0 is used.
  if (param->threads == 0) {
    Logging::Warn("--threads resolved to 0 (hardware_concurrency unavailable). Defaulting to 1.");
    param->threads = 1;
  }
  // Clamp worker_threads_per_input to the thread count so the default (4) works with any --threads value.
  if (param->worker_threads_per_input > param->threads) {
    Logging::Warn("--workers-per-reads-file ({}) exceeds --threads ({}). Clamping to {}.",
                  param->worker_threads_per_input,
                  param->threads,
                  param->threads);
    param->worker_threads_per_input = param->threads;
  }
  // Make sure that min_trimmed_read_len less than or equal to min_read_len.
  // It is impossible trimmed reads to be longer than raw reads.
  if (param->min_trimmed_read_len > param->min_read_len) {
    Logging::Warn(
        "--min-trimmed-read-length ({}) is greater than --min-read-length ({}). "
        "Trimmed reads cannot be longer than raw reads so we will set --min-trimmed-read-length to "
        "--min-read-length.",
        param->min_trimmed_read_len,
        param->min_read_len);
    param->min_trimmed_read_len = param->min_read_len;
  }

  // if the output directory is a file, then throw an error
  if (fs::exists(param->out_dir) && !fs::is_directory(param->out_dir)) {
    throw CLI::ValidationError(fmt::format("Output path ({}) is not a directory.", param->out_dir));
  }

  // create directory if it does not exist
  if (!fs::exists(param->out_dir)) {
    fs::create_directories(param->out_dir);
  }
  // if it exists, check if we can write to it
  if (fs::exists(param->out_dir)) {
    // throw an error if unable to write to the directory
    file::CheckFileIsWritable(param->out_dir);
  } else {
    // not created so throwing an error
    throw CLI::ValidationError(
        fmt::format("Output directory ({}) does not exist and could not be created.", param->out_dir));
  }
}

static void PreMainProcessing(const CLI::App* const app,
                              const std::shared_ptr<DemuxAndTrimParam>& param,
                              const DeprecatedDemuxCliOptionsPtr& deprecated) {
  ApplyDeprecatedOptions(app, param, deprecated);
  ValidateOptions(param);
  param->command_line_info = cli::GetCommandLineInfo(app);
}

DemuxCliCallbacks CreateCliCallbacks() {
  const auto deprecated = std::make_shared<DeprecatedDemuxCliOptions>();
  return {
      .define_options =
          [deprecated](CLI::App* const app, const std::shared_ptr<DemuxAndTrimParam>& param) {
            DefineGeneralOptions(app, param, deprecated);
            DefineSimplexOptions(app, param);
            DefineDuplexOptions(app, param);
          },
      .pre_callback =
          [deprecated](const CLI::App* const app, const std::shared_ptr<DemuxAndTrimParam>& param) {
            PreMainProcessing(app, param, deprecated);
          },
  };
}

}  // namespace xoos::demux
