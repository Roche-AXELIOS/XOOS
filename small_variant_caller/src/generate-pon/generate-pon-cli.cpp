#include "generate-pon-cli.h"

#include <string>

#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/file-permission-validator.h>
#include <xoos/cli/validators/nonempty-file-validator.h>

#include "core/cli-option-names.h"
#include "util/cli-util.h"

namespace xoos::svc {

constexpr auto* kDefaultOutputFile = "pon.vcf.gz";
constexpr f32 kDefaultMinDuplexAf = 0.01F;
constexpr f32 kDefaultMaxDuplexAf = 0.1F;
constexpr f32 kDefaultMinMapqMean = 10.0F;

namespace {
/**
 * @brief Extract a list of file paths from a text file.
 * @param argument_name CLI argument name (for error messages).
 * @param list_file_name Path to the text file containing one file path per line.
 * @param files Output vector populated with the parsed file paths.
 * @throws cli::ValidationError if any listed file is not readable.
 */
void ParseFileList(const std::string& argument_name, const fs::path& list_file_name, vec<fs::path>& files) {
  const auto results = cli::CheckFileListReadable(list_file_name);
  if (std::holds_alternative<std::string>(results)) {
    throw cli::ValidationError("Invalid file list for '{}'; {}", argument_name, std::get<std::string>(results));
  }
  files = std::get<vec<fs::path>>(results);
}
}  // namespace

void DefineOptionsGeneratePon(CLI::App* const app, const GeneratePonParamPtr& params) {
  auto* const feature_list_opt = app->add_option_function<fs::path>(
      cli_opt_name::kFeatureList,
      [&params](const fs::path& value) { ParseFileList(cli_opt_name::kFeatureList, value, params->feature_files); },
      "Text file containing a list of BAM feature file paths, one per line. "
      "Each file is expected to be a TSV output from compute_bam_features with the pon workflow.");
  feature_list_opt->required();
  feature_list_opt->check(cli::NonEmptyFileValidator());

  app->add_option(cli_opt_name::kOutputFile, params->output_file, "Output PON VCF file path")
      ->default_val(kDefaultOutputFile)
      ->check(CLI::NonexistentPath);

  app->add_option(cli_opt_name::kMinDuplexAf,
                  params->min_duplex_af,
                  "Minimum duplex_af to include a variant (variants below this are excluded)")
      ->default_val(kDefaultMinDuplexAf)
      ->check(kCliRangeFraction);

  app->add_option(cli_opt_name::kMaxDuplexAf,
                  params->max_duplex_af,
                  "Maximum duplex_af to include a variant (variants above this are excluded)")
      ->default_val(kDefaultMaxDuplexAf)
      ->check(kCliRangeFraction);

  app->add_option(cli_opt_name::kMinMapqMean,
                  params->min_mapq_mean,
                  "Minimum mapq_mean to include a variant (variants below this are excluded)")
      ->default_val(kDefaultMinMapqMean)
      ->check(CLI::NonNegativeNumber);

  app->callback([&params]() {
    if (params->min_duplex_af > params->max_duplex_af) {
      throw cli::ValidationError("{} ({}) must be <= {} ({})",
                                 cli_opt_name::kMinDuplexAf,
                                 params->min_duplex_af,
                                 cli_opt_name::kMaxDuplexAf,
                                 params->max_duplex_af);
    }
  });

  cli::AddThreadCountOption(app, cli_opt_name::kThreads, params->threads);
  AddWarnAsErrorOption(app);
}

}  // namespace xoos::svc

namespace xoos::svc::generate_pon {

void PreCallback(const CLI::App* const app, const GeneratePonParamPtr& params) {
  params->command_line = cli::GetCommandLineInfo(app);
}

}  // namespace xoos::svc::generate_pon
