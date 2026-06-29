#pragma once

#include <memory>

#include <xoos/cli/cli.h>
#include <xoos/cli/enum-option-util.h>

#include "compute-bam-features.h"
#include "core/cli-option-names.h"
#include "util/cli-util.h"

namespace xoos::svc {
using ComputeBamFeaturesCliParamsPtr = std::shared_ptr<ComputeBamFeaturesCliParams>;
}  // namespace xoos::svc

namespace xoos::svc::compute_bam_features {
/**
 * @brief Add shared CLI options for a subcommand with workflow-specific default values.
 * @details This template function is intended to be used in subcommands of compute_bam_features and other SVC
 * submodules that require BAM feature extraction, such as filter_variants. Vector of added CLI options is returned for
 * further processing, if needed.
 * @tparam Params Type of CLI parameters
 * @param sub Subcommand application pointer where CLI options are to be added
 * @param params CLI parameters shared pointer to store option values
 * @param defaults Config containing workflow-specific default values for the options
 * @return Vector of all CLI options added to the subcommand
 */
template <typename Params>
vec<CLI::Option*> AddSharedOptions(CLI::App* const sub, std::shared_ptr<Params>& params, const SVCConfig& defaults) {
  vec<CLI::Option*> options;

  // Add options with workflow-specific default values

  options.push_back(
      sub->add_option(cli_opt_name::kMaxVariantsPerRead,
                      params->max_read_variant_count,
                      "Max number of variants allowed per read (inclusive); `0` can also turn off this option")
          ->default_val(defaults.max_variants_per_read)
          ->check(CLI::NonNegativeNumber));

  options.push_back(
      sub->add_option(
             cli_opt_name::kMaxVariantsPerReadNormalized,
             params->max_read_variant_count_normalized,
             "Max number of variants allowed per read, normalized by alignment length (inclusive); `0` can also "
             "turn off this option")
          ->default_val(defaults.max_variants_per_read_normalized)
          ->check(kCliRangeFraction));

  options.push_back(sub->add_option(
      cli_opt_name::kSkipVariantsVcf,
      params->skip_variants_vcf,
      "VCF containing variants not counted by `--max-variants-per-read` or `--max-variants-per-read-normalized`"));
  CheckIndexedVcfFile(options.back());

  options.push_back(cli::AddEnumOption(sub,
                                       cli_opt_name::kSequencingProtocol,
                                       params->sequencing_protocol,
                                       "the sequencing protocol used to generate the input data",
                                       defaults.sequencing_protocol));

  options.push_back(cli::AddEnumOption(sub,
                                       cli_opt_name::kDecodeYc,
                                       params->decode_yc,
                                       "decode YC tags within input BAM file(s) using the specified method",
                                       defaults.decode_yc));

  options.push_back(cli::AddEnumOption(sub,
                                       cli_opt_name::kMinBaseType,
                                       params->min_base_type,
                                       "minimum base type in duplex reads for variant support",
                                       defaults.min_base_type));

  options.push_back(
      cli::AddEnumOption(sub,
                         cli_opt_name::kFilterHomopolymer,
                         params->filter_homopolymer,
                         "skip variant adjacent to homopolymer that spans beyond the alignment's end position",
                         defaults.filter_homopolymer));

  // Base quality and mapping quality are stored as 8-bit unsigned int in the CLI params struct and the config struct.
  // Static-cast the default value to a 16-bit unsigned int to bypass the conversion error in the CLI11 library, where
  // the 8-bit value being treated as a char instead of an integer.
  options.push_back(sub->add_option(cli_opt_name::kMinBq,
                                    params->min_bq,
                                    "Minimum alignment base quality required to support a variant (inclusive)")
                        ->default_val(static_cast<u16>(defaults.min_bq))
                        ->check(kCliRangeBaseq));

  options.push_back(sub->add_option(cli_opt_name::kMinMapq,
                                    params->min_mapq,
                                    "Minimum alignment mapping quality required to support a variant (inclusive)")
                        ->default_val(static_cast<u16>(defaults.min_mapq))
                        ->check(kCliRangeMapq));

  options.push_back(sub->add_option(cli_opt_name::kMinDist,
                                    params->min_allowed_distance_from_end,
                                    "Minimum distance of variant from fragment alignment end (inclusive)")
                        ->default_val(defaults.min_dist)
                        ->check(CLI::NonNegativeNumber));

  options.push_back(sub->add_option(cli_opt_name::kMinFamilySize,
                                    params->min_family_size,
                                    "Minimum cluster size required for an alignment to support a variant (inclusive)")
                        ->default_val(defaults.min_family_size)
                        ->check(CLI::NonNegativeNumber));

  options.push_back(sub->add_option(cli_opt_name::kMinHomopolymerLength,
                                    params->min_homopolymer_length,
                                    "Minimum length of homopolymer in reference")
                        ->default_val(defaults.min_homopolymer_length)
                        ->check(CLI::NonNegativeNumber));

  return options;
}

/**
 * @brief Apply workflow config to CLI parameters if they were not specified via the subcommand CLI.
 * @details This template function is intended to be used in the pre-callback functions of compute_bam_features and
 * other SVC submodules that require BAM feature extraction, such as filter_variants.
 * @tparam Params Type of CLI parameters
 * @param sub Subcommand CLI application pointer to check which options were specified via CLI
 * @param params CLI parameters shared pointer to apply defaults to
 */
template <typename Params>
void ApplyConfig(const cli::ConstAppPtr sub, const std::shared_ptr<Params>& params) {
  // Only set params fields with workflow-specific config preset if they were not specified via CLI
  if (sub->count(cli_opt_name::kMinMapq) == 0) {
    params->min_mapq = params->config.min_mapq;
  }
  if (sub->count(cli_opt_name::kMinBq) == 0) {
    params->min_bq = params->config.min_bq;
  }
  if (sub->count(cli_opt_name::kMinDist) == 0) {
    params->min_allowed_distance_from_end = params->config.min_dist;
  }
  if (sub->count(cli_opt_name::kMinFamilySize) == 0) {
    params->min_family_size = params->config.min_family_size;
  }
  if (sub->count(cli_opt_name::kFilterHomopolymer) == 0) {
    params->filter_homopolymer = params->config.filter_homopolymer;
  }
  if (sub->count(cli_opt_name::kMinHomopolymerLength) == 0) {
    params->min_homopolymer_length = params->config.min_homopolymer_length;
  }
  if (sub->count(cli_opt_name::kSequencingProtocol) == 0) {
    params->sequencing_protocol = params->config.sequencing_protocol;
  }
  if (sub->count(cli_opt_name::kDecodeYc) == 0) {
    params->decode_yc = params->config.decode_yc;
  }
  if (sub->count(cli_opt_name::kMinBaseType) == 0) {
    params->min_base_type = params->config.min_base_type;
  }
}

/**
 * @brief Define CLI options for the main application.
 * @param app Main application pointer where CLI options will be defined.
 * @param params Shared pointer to CLI parameters to store parsed option values
 */
void DefineOptions(CLI::App* app, ComputeBamFeaturesCliParamsPtr& params);

/**
 * @brief Preprocess the CLI options and validate the parameters before running the main application.
 * @param app Main application pointer where the callback is triggered
 * @param params Shared pointer to CLI parameters to be validated
 * @throws CLI::ValidationError if the parameters are invalid
 */
void PreCallback(cli::ConstAppPtr app, const ComputeBamFeaturesCliParamsPtr& params);
}  // namespace xoos::svc::compute_bam_features
