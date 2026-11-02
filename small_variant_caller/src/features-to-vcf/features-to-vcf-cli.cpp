#include "features-to-vcf-cli.h"

#include <string>

#include <xoos/cli/validators/indexed-fasta-validator.h>
#include <xoos/cli/validators/nonempty-file-validator.h>

#include "core/cli-option-names.h"
#include "util/cli-util.h"

namespace xoos::svc {

namespace {
const std::string kDefaultOutputFile = "features.vcf.gz";
const std::string kDefaultSampleName = "SAMPLE";
}  // namespace

void DefineOptionsFeaturesToVcf(CLI::App* const app, const FeaturesToVcfParamPtr& params) {
  auto* const input_opt =
      app->add_option(cli_opt_name::kInputFile, params->input_file, "Path to input BAM features TSV file")->required();
  input_opt->check(cli::NonEmptyFileValidator());

  auto* const genome_opt =
      app->add_option(cli_opt_name::kGenome, params->genome, "Path to indexed FASTA reference genome")->required();
  genome_opt->check(cli::IndexedFastaFileValidator());

  app->add_option(cli_opt_name::kOutputFile, params->output_file, "Path to output VCF file")
      ->default_val(kDefaultOutputFile);
  app->add_option(cli_opt_name::kSampleName, params->sample, "Sample name for the VCF header")
      ->default_val(kDefaultSampleName);
  AddWarnAsErrorOption(app);
}

}  // namespace xoos::svc
