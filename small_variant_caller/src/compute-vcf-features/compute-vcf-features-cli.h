#pragma once

#include <memory>

#include <xoos/cli/cli.h>

#include "compute-vcf-features.h"
#include "core/cli-option-names.h"
#include "util/cli-util.h"

namespace xoos::svc {
using ComputeVcfFeaturesParamPtr = std::shared_ptr<ComputeVcfFeaturesParam>;
}  // namespace xoos::svc

namespace xoos::svc::compute_vcf_features {
/**
 * @brief Add shared CLI options for a subcommand with workflow-specific default values.
 * @details This template function is intended to be used in subcommands of compute_vcf_features and other SVC
 * submodules that require VCF feature extraction, such as filter_variants. Vector of added CLI options is returned for
 * further processing, if needed.
 * @tparam Params Type of CLI parameters
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params CLI parameters shared pointer to store option values
 * @return Vector of all CLI options added
 */
template <typename Params>
vec<CLI::Option*> AddSharedOptions(CLI::App* const sub, std::shared_ptr<Params>& params) {
  vec<CLI::Option*> options;

  options.push_back(sub->add_option_function<fs::path>(
      cli_opt_name::kInterestRegions,
      [&params](const fs::path& value) { params->interest_regions = GetChromIntervalMap(value); },
      "Path to a BED file of regions for VCF feature `at_interest_region`"));
  CheckBedFile(options.back());

  options.push_back(
      sub->add_option(cli_opt_name::kPopAfVcf,
                      params->pop_af_vcf,
                      "Path to GATK gnomAD population allele frequency VCF (i.e. af-only-gnomad.hg38.vcf.gz)"));
  CheckIndexedVcfFile(options.back());

  return options;
}

/**
 * @brief Define CLI options for the main application.
 * @param app Main application pointer where CLI options will be defined.
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void DefineOptions(CLI::App* app, ComputeVcfFeaturesParamPtr& params);

/**
 * @brief Preprocess the CLI options and validate the parameters before running the main application.
 * @param app Main application pointer where the callback is triggered
 * @param params Shared pointer to CLI parameters to be validated
 * @throws CLI::ValidationError if the parameters are invalid
 */
void PreCallback(cli::ConstAppPtr app, const ComputeVcfFeaturesParamPtr& params);
}  // namespace xoos::svc::compute_vcf_features
