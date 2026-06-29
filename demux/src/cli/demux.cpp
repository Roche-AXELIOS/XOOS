#include "demux.h"

#include <xoos/cli/enum-option-util.h>
#include <xoos/cli/file-option-util.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/enum/enum-util.h>
#include <xoos/io/metadata-util.h>
#include <xoos/util/file-functions.h>

#include <fstream>
#include <string>
#include <vector>

#include "io/sequence-reader.h"
#include "task/batch.h"

namespace xoos::demux {

using cli::AddEnumOption;
using cli::AddOutputFileOption;
using cli::AddThreadCountOption;

/**
 * Return all sequence files from a directory recursive manner
 */
static std::vector<fs::path> RecursiveFindSequenceFiles(const fs::path& input) {
  std::vector<fs::path> stack{input};

  std::vector<fs::path> files;
  while (!stack.empty()) {
    auto current = stack.back();
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

std::vector<fs::path> ReadInputFileList(const fs::path& input_file_list) {
  std::ifstream ifs{input_file_list};
  std::vector<fs::path> input_files;
  std::string line;
  while (std::getline(ifs, line)) {
    if (!fs::exists(line)) {
      throw CLI::ValidationError(fmt::format("Input does not exist: '{}'", line));
    }
    input_files.emplace_back(line);
  }
  return input_files;
}

std::vector<fs::path> ExpandInputFileList(const std::vector<fs::path>& input_files) {
  std::vector<fs::path> expanded_input_files;
  for (const auto& input_file : input_files) {
    auto files{RecursiveFindSequenceFiles(input_file)};
    expanded_input_files.insert(expanded_input_files.end(), files.cbegin(), files.cend());
  }
  return expanded_input_files;
}

void DefineOptions(cli::AppPtr app, std::shared_ptr<DemuxAndTrimParam>& param) {
  app->add_option_function<std::vector<fs::path>>(
         "-i,--input",
         [&param](const std::vector<fs::path>& inputs) -> void {
           const auto expanded{ExpandInputFileList(inputs)};
           std::copy(std::cbegin(expanded), std::cend(expanded), std::back_inserter(param->inputs));
         },
         "Inputs containing untrimmed sequences in either FASTQ, FASTQ.gz, or RDB format. "
         "May specify a directory (reading all relevant files in directory). May be specified multiple times, or as a "
         "space separated list.")
      ->check(CLI::ExistingPath);
  app->add_option_function<fs::path>(
         "--input-files-list",
         [&param](const std::string& input_file_list) -> void {
           auto inputs = ExpandInputFileList(ReadInputFileList(input_file_list));
           std::copy(std::cbegin(inputs), std::cend(inputs), std::back_inserter(param->inputs));
         },
         "File listing the input file paths instead of setting them via -i/--input option.")
      ->check(CLI::ExistingFile);
  app->add_option("-o,--out-dir", param->out_dir, "Output directory")->default_val(param->out_dir);
  app->add_flag("--overwrite", param->overwrite, "Overwrites existing non-empty output directory and its contents.")
      ->default_val(param->overwrite);
  app->add_option("-b,--batch-size", param->batch_size, "Number of reads to process in each batch.")
      ->default_val(param->batch_size);
  AddThreadCountOption(app, "--threads", param->threads);
  AddEnumOption(app, "--read-length-mode", param->read_length_mode,
                "Read filtering type for demultiplexing simplex reads.", param->read_length_mode);
  app->add_option("-p,--adapter-design-bundle", param->adapter_design_bundle,
                  "Path to ZIP file or directory containing the adapter architecture, search strategy, and default "
                  "adapter design")
      ->default_val(param->adapter_design_bundle)
      ->check(CLI::ExistingPath);
  app->add_option("-n,--adapter-design-name", param->adapter_design_name,
                  "The adapter design name to be loaded from the bundle. Defaults to SBX-D.")
      ->check(CLI::Validator(
          // Validation mostly to control methylation mode for now
          [&param](const std::string& input) -> std::string {
            if (input.empty()) {
              return "adapter-design-name must not be empty";
            }
            if (input == "SBX-DM") {
              param->methylation = true;
            }
            return {};
          },
          "valid-adapter-design-name"))
      ->default_val(param->adapter_design_name);
  app->add_option(
         "--sample-sheet", param->sample_sheet,
         "If specified only look for the adapter SIDs in sample sheet. Otherwise adapter design bundle sheet is used.")
      ->check(CLI::ExistingFile);
  app->add_option(
         "--compression-level", param->compression_level,
         "Compression level for output FASTQ if using gzip (1-9) or zstd (1-19) compression. Higher level increases "
         "compression at the cost of speed.  If compression type is set to none, this option is ignored and "
         "output will be uncompressed.  Will choose default for both algorithms if not specified.")
      ->check(CLI::Range(kMinCompressionLevel, kMaxZstdCompressionLevel));
  app->add_option("--writing-threads-per-sample", param->writing_threads_per_sample,
                  "Max number of writer threads to be used for demux output (1-8=number of threads)")
      ->check(CLI::Range(size_t{1}, BatchTask::kMaxNumberOfSinkWorkers))
      ->default_val(param->writing_threads_per_sample);
  AddEnumOption(
      app, "--compression-type", param->compression_type,
      "Compression type for output FASTQs. If set to none, output will ignore the compression level parameter.",
      param->compression_type);
  app->add_option("--length-distribution-report-max", param->length_distribution_report_max,
                  "The maximum size the largest bucket in the length histogram distribution report can be.")
      ->default_val(param->length_distribution_report_max);
  app->add_option("--min-read-len", param->min_read_len, "Filter raw reads shorter than this length.")
      ->check(CLI::PositiveNumber)
      ->default_val(param->min_read_len);
  app->add_option("--min-trimmed-read-len", param->min_trimmed_read_len,
                  "Filter trimmed reads shorter than this length.")
      ->check(CLI::PositiveNumber)
      ->default_val(param->min_trimmed_read_len);
  app->add_option("--max-error-rate-percent", param->max_error_rate_percent,
                  "Filter duplex reads with a consensus region with error rate larger than this percentage.")
      ->default_val(param->max_error_rate_percent);
  app->add_option("--stop-after-min-concordant-duplex-bases", param->min_concord_dp_bases,
                  "Minimum concordant duplex bases per sample for early stopping.");
  app->add_option("--strand-critical-phred", param->strand_critical,
                  "In duplex strand detection mode, this is the Phred scaled FPR critical value "
                  "-10*log10(critical_value). Low values will classify faster at cost of accuracy.")
      ->default_val(param->strand_critical);
  app->add_option_function<fs::path>(
      "--strand-detect",
      [&param](const fs::path& filename) -> void {
        param->strand_detector.emplace(filename);
        auto error_rates = param->strand_detector.value().FalsePositiveRates();
        Logging::Info("Strand Detection Mode. Filter FPRs fw: {:.5f} rv: {:.5f}", error_rates.first,
                      error_rates.second);
        double average_error_rate = (error_rates.first + error_rates.second) / 2.0;
        // if the strand detector is specified, we need to set the strand classifier
        double rescaled_critical_value = std::pow(10, -param->strand_critical / 10.0);
        param->strand_classifier.emplace(param->strand_detector.value(), average_error_rate, rescaled_critical_value);
      },
      "Enable genome strand detection with the specified file created by demux-strand-index.");
  app->add_flag("--output-failed-reads", param->output_failed_reads,
                "Output unaltered failed reads into files with prefix 'raw_failed'.")
      ->default_val(param->output_failed_reads);
}

void ValidateOptions(const std::shared_ptr<DemuxAndTrimParam>& param) {
  // Make sure that we have at least one input file.
  if (param->inputs.empty()) {
    throw CLI::ValidationError("Must provide at least one input");
  }
  // Make sure that min_trimmed_read_len less than or equal to min_read_len.
  // It is impossible trimmed reads to be longer than raw reads.
  if (param->min_trimmed_read_len > param->min_read_len) {
    Logging::Warn(
        "--min-trimmed-read-len ({}) is greater than --min-read-len ({}). "
        "Trimmed reads cannot be longer than raw reads so we will set --min-trimmed-read-len to --min-read-len.",
        param->min_trimmed_read_len, param->min_read_len);
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

void PreMainProcessing(cli::ConstAppPtr app, const std::shared_ptr<DemuxAndTrimParam>& param,
                       const std::string& version) {
  ValidateOptions(param);
  AddVersionAndCommandLineComment(app, param, version);
}

cli::PreCallback<DemuxAndTrimParam> CreatePreCallback(const std::string& version) {
  return [&version](cli::ConstAppPtr app, const std::shared_ptr<DemuxAndTrimParam>& param) -> void {
    PreMainProcessing(app, param, version);
  };
}

// This adds version and commandline to the comments field of the param
// This handles parent or subcommands appropriately if available.
void AddVersionAndCommandLineComment(const cli::ConstAppPtr& app, const std::shared_ptr<DemuxAndTrimParam>& param,
                                     const std::string& version) {
  const std::string cli_args = cli::RenderFullCli(app);
  io::AddVersionAndCommandLineComment(param->comments, version, cli_args);
}

}  // namespace xoos::demux
