#include "detect-sample-type-cli.h"

#include <fmt/core.h>

#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/bed-validator.h>
#include <xoos/cli/validators/nonempty-file-validator.h>
#include <xoos/enum/enum-util.h>
#include <xoos/log/logging.h>

#include "core/cli-option-names.h"
#include "util/cli-util.h"

namespace xoos::svc {

void DefineOptionsDetectSampleType(CLI::App* const app, const DetectSampleTypeParamPtr& params) {
  auto* const vcf_opt =
      app->add_option(cli_opt_name::kVcfInput, params->vcf_file, "Path to input tumor-normal VCF file")->required();
  vcf_opt->check(cli::NonEmptyFileValidator());

  auto* const target_bed_opt = app->add_option(
      cli_opt_name::kTargetBed,
      params->target_bed,
      "Optional uncompressed or gzip/BGZF-compressed (.gz) BED file limiting sample-type detection to target regions");
  target_bed_opt->check(cli::BedFileValidator());

  app->add_option(cli_opt_name::kMaxAf,
                  params->options.max_af,
                  "Maximum allele frequency for SNV error observations used by sample-type detection")
      ->default_val(kDefaultMaxAfForSampleTypeDetection)
      ->check(kCliRangeFraction);

  app->add_option(cli_opt_name::kMaxFfpePhred,
                  params->options.max_ffpe_phred,
                  "Maximum VCF phred for classifying sample type as FFPE")
      ->default_val(kDefaultMaxFfpePhredForSampleTypeDetection)
      ->check(CLI::NonNegativeNumber);

  cli::AddThreadCountOption(app, cli_opt_name::kThreads, params->options.threads);
  AddWarnAsErrorOption(app);
}

void DetectSampleType(const DetectSampleTypeParam& params) {
  const auto metrics = ComputeSampleTypeMetrics(params.vcf_file, params.target_bed, params.options);
  const auto sample_type = ClassifySampleType(metrics.phred, params.options.max_ffpe_phred);
  Logging::Info(
      "Inferred sample type '{}' from VCF phred {:.4f}, max FFPE phred {}, median DP {}, error count {}, "
      "total length {}",
      enum_util::FormatEnumName(sample_type),
      metrics.phred,
      params.options.max_ffpe_phred,
      metrics.median_dp,
      metrics.error_count,
      metrics.total_length);
  fmt::print("{}\n", enum_util::FormatEnumName(sample_type));
}

}  // namespace xoos::svc
