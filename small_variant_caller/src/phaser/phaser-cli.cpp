#include "phaser-cli.h"

#include <array>

#include <fmt/format.h>

#include <xoos/cli/cli.h>
#include <xoos/cli/thread-count-option-util.h>
#include <xoos/cli/validators/indexed-bam-validator.h>
#include <xoos/cli/validators/indexed-fasta-validator.h>
#include <xoos/cli/validators/nonempty-file-validator.h>
#include <xoos/enum/enum-util.h>

#include "CLI/CLI.hpp"
#include "core/cli-option-names.h"
#include "core/workflow.h"
#include "util/cli-util.h"

namespace xoos::svc::phaser {

using cli::AddThreadCountOption;
using enum Workflow;

namespace {

constexpr auto* kDefaultPhasedVcfOutput = "output.phased.vcf.gz";

void AddCommonOptions(CLI::App* const sub, const PhaserParamPtr& params, const bool require_tumor_sample_name) {
  auto* const bam_opt =
      sub->add_option(cli_opt_name::kBamInput, params->bam_path, "Path to Mutect2 bamout BAM file")->required();
  bam_opt->check(cli::IndexedBamFileValidator());

  auto* const vcf_opt =
      sub->add_option(cli_opt_name::kVcfInput, params->vcf_path, "Path to XOOS somatic VCF")->required();
  vcf_opt->check(cli::NonEmptyFileValidator());

  sub->add_option(cli_opt_name::kOutputDir, params->output_dir, "Path to output directory")->default_val(".");
  sub->add_option(cli_opt_name::kVcfOutput,
                  params->vcf_output,
                  "Output VCF filename, relative to --output-dir (default: output.phased.vcf.gz)")
      ->default_val(kDefaultPhasedVcfOutput);

  auto* const ref_opt = sub->add_option(
      cli_opt_name::kGenome,
      params->ref_path,
      "Optional reference FASTA used as the backbone when merging phased variants; the input BAM is used as a "
      "fallback");
  ref_opt->check(cli::IndexedFastaFileValidator());

  if (require_tumor_sample_name) {
    sub->add_option(cli_opt_name::kTumorSampleName,
                    params->tumor_rg,
                    "Tumor read group ID in BAM for the tumor-normal-wgs workflow")
        ->required();
  }

  sub->add_option(cli_opt_name::kMaxPhaseDistance,
                  params->max_phase_distance,
                  "Maximum distance between variants to consider for phasing")
      ->default_val(params->max_phase_distance)
      ->check(CLI::NonNegativeNumber);
  sub->add_option(
         cli_opt_name::kMinSharedReads, params->min_shared_reads, "Minimum reads spanning both variant positions")
      ->default_val(params->min_shared_reads)
      ->check(CLI::NonNegativeNumber);
  sub->add_option(cli_opt_name::kMinCisReads, params->min_cis_reads, "Minimum reads supporting cis configuration")
      ->default_val(params->min_cis_reads)
      ->check(CLI::NonNegativeNumber);
  sub->add_option(cli_opt_name::kMinCisFraction,
                  params->min_cis_fraction,
                  "Minimum fraction of cis reads among informative reads")
      ->default_val(params->min_cis_fraction)
      ->check(kCliRangeFraction);
  sub->add_flag(cli_opt_name::kReportMnv,
                params->report_mnv,
                "Annotate phased variants with an MNVTAG INFO field (use --merge-mnv/--merge-complex "
                "to emit merged MNV records)");
  sub->add_flag(cli_opt_name::kMergeMnv, params->merge_mnv, "Merge phased adjacent SNVs into single MNV records");
  sub->add_flag(cli_opt_name::kMergeComplex,
                params->merge_complex,
                "Merge all phased variants (including indels) into single records");
  sub->add_flag(cli_opt_name::kKeepComponents,
                params->keep_components,
                "When merging, keep original component variant lines alongside merged records");
  sub->add_flag(
      cli_opt_name::kIncludeFiltered, params->include_filtered, "Include filtered variants in phasing analysis");
  sub->add_flag(cli_opt_name::kRescueFiltered,
                params->rescue_filtered,
                "Phase PASS+non-PASS pairs; rescue non-PASS variants if they phase with a PASS variant");
  sub->add_flag(cli_opt_name::kSkipInsSnv,
                params->skip_ins_snv,
                "Skip phasing of insertion+SNV pairs (reduces false positives near insertions)");
  sub->add_option(cli_opt_name::kMinBaseQualityOpt,
                  params->min_base_quality,
                  "Minimum base quality at variant position to count a read")
      ->default_val(params->min_base_quality)
      ->check(kCliRangeBaseq);
  AddThreadCountOption(sub, cli_opt_name::kThreads, params->threads);
  AddWarnAsErrorOption(sub);
}

}  // namespace

void DefineOptions(CLI::App* const app, const PhaserParamPtr& params) {
  constexpr std::array kSupportedWorkflows = {kTumorNormalWgs, kTumorOnlyWgs};
  for (const Workflow workflow : kSupportedWorkflows) {
    const std::string name = enum_util::FormatEnumName(workflow);
    const std::string description = fmt::format("Phase variants for the {} workflow", name);
    auto* const sub = app->add_subcommand(name, description)->fallthrough();
    AddCommonOptions(sub, params, workflow == kTumorNormalWgs);
  }
  app->require_subcommand(1, 1);
}

void PreCallback(const CLI::App* const app, const PhaserParamPtr& params) {
  if (app->got_subcommand(enum_util::FormatEnumName(kTumorNormalWgs))) {
    params->workflow = kTumorNormalWgs;
  } else if (app->got_subcommand(enum_util::FormatEnumName(kTumorOnlyWgs))) {
    params->workflow = kTumorOnlyWgs;
  }
  params->output_path = params->vcf_output.is_absolute() ? params->vcf_output : params->output_dir / params->vcf_output;
  if (fs::exists(params->output_path)) {
    throw CLI::ValidationError(fmt::format("Output VCF file '{}' already exists", params->output_path.string()));
  }
  if (params->keep_components && !(params->merge_mnv || params->merge_complex)) {
    throw CLI::ValidationError(cli_opt_name::kKeepComponents, "requires --merge-mnv or --merge-complex to be enabled");
  }
  params->command_line = cli::GetCommandLineInfo(app);
}

}  // namespace xoos::svc::phaser
