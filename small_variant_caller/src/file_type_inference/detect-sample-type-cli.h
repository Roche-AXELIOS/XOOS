#pragma once

#include <memory>
#include <optional>

#include <xoos/cli/cli.h>
#include <xoos/types/fs.h>

#include "file_type_inference/detect-sample-type.h"

namespace xoos::svc {

/**
 * @brief Parsed parameters for the detect_sample_type CLI.
 * @details The CLI stores values directly in this struct, then DetectSampleType() passes the VCF, optional
 * BED, and detection options to the shared sample-type detection implementation.
 */
struct DetectSampleTypeParam {
  /// Input VCF or VCF.GZ used to compute sample-type metrics.
  fs::path vcf_file;
  /// Optional BED target file restricting both evaluated records and denominator length.
  std::optional<fs::path> target_bed;
  /// Detection thresholds and threading configuration.
  DetectSampleTypeOptions options;
};

using DetectSampleTypeParamPtr = std::shared_ptr<DetectSampleTypeParam>;

/**
 * @brief Register CLI options for the detect_sample_type subcommand.
 * @details Adds required VCF input, optional BED targeting, detection threshold, threading, and warning
 * handling options. Values are written directly into @p params for use by DetectSampleType().
 * @param app CLI11 application/subcommand receiving the detection options.
 * @param params Shared parameter object populated by CLI parsing.
 */
void DefineOptionsDetectSampleType(CLI::App* app, const DetectSampleTypeParamPtr& params);

/**
 * @brief Execute sample-type detection from parsed CLI parameters.
 * @details Computes the sample-type metrics, applies the FFPE phred threshold, logs the supporting metrics,
 * and writes the inferred sample type to stdout for command-line use.
 * @param params Parsed detect_sample_type parameters.
 */
void DetectSampleType(const DetectSampleTypeParam& params);

}  // namespace xoos::svc
