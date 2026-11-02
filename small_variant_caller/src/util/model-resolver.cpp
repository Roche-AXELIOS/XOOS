#include "model-resolver.h"

#include <fmt/format.h>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include "core/cli-option-names.h"

namespace xoos::svc {

using enum AlignerType;

namespace {
/**
 * @brief Convert an aligner enum to the model filename token.
 * @details Rejects AlignerType::kCustom because custom aligners require explicit model paths and should
 * not be resolved through bundled model naming conventions.
 * @param aligner Aligner type to render.
 * @return Formatted aligner name used in bundled model filenames.
 * @throws xoos::error::Error when @p aligner is AlignerType::kCustom.
 */
std::string AlignerToModelString(const AlignerType aligner) {
  if (aligner == kCustom) {
    throw error::Error(fmt::format("Cannot resolve model for {} aligner — provide {} explicitly",
                                   enum_util::FormatEnumName(aligner),
                                   cli_opt_name::kModel));
  }
  return enum_util::FormatEnumName(aligner);
}

/**
 * @brief Format an optional model-threshold value for logging.
 * @details Keeps log messages explicit by rendering absent thresholds as "unset" rather than an empty value.
 * @param threshold Optional threshold value.
 * @return String representation of the threshold or "unset".
 */
std::string FormatOptionalThreshold(const std::optional<f32>& threshold) {
  return threshold.has_value() ? fmt::format("{}", *threshold) : "unset";
}
}  // namespace

fs::path ResolveTumorNormalModel(const SampleType sample_type, const AlignerType aligner) {
  if (sample_type == SampleType::kAutoDetect) {
    throw error::Error("Cannot resolve tumor-normal model while sample type is '{}'",
                       enum_util::FormatEnumName(SampleType::kAutoDetect));
  }
  // No pre-trained giraffe models exist for tumor-normal-wgs — reject rather than silently
  // resolving a mismatched model.
  if (aligner == kGiraffe) {
    throw error::Error(
        fmt::format("Aligner '{}' is not supported for the tumor-normal-wgs workflow because no pre-trained {} model "
                    "is available. Use '{} {}', or '{} {}' with an explicit {}.",
                    enum_util::FormatEnumName(kGiraffe),
                    enum_util::FormatEnumName(kGiraffe),
                    cli_opt_name::kAligner,
                    enum_util::FormatEnumName(kBwa),
                    cli_opt_name::kAligner,
                    enum_util::FormatEnumName(AlignerType::kCustom),
                    cli_opt_name::kModel));
  }
  const auto sample_type_str = enum_util::FormatEnumName(sample_type);
  const auto path =
      fmt::format("/resources/model-tumor-normal-wgs-{}-{}.txt.gz", sample_type_str, AlignerToModelString(aligner));

  Logging::Info("Resolved model: {} based on aligner {} and sample_type {}",
                path,
                enum_util::FormatEnumName(aligner),
                sample_type_str);
  return fs::path{path};
}

fs::path ResolveTumorOnlyWgsModel(const SampleType sample_type, const AlignerType aligner) {
  using enum SampleType;
  // No pre-trained giraffe models exist for tumor-only-wgs — reject rather than silently
  // resolving a mismatched model.
  if (aligner == kGiraffe) {
    throw error::Error(
        fmt::format("Aligner '{}' is not supported for the tumor-only-wgs workflow because no pre-trained {} model "
                    "is available. Use '{} {}', or '{} {}' with an explicit {}.",
                    enum_util::FormatEnumName(kGiraffe),
                    enum_util::FormatEnumName(kGiraffe),
                    cli_opt_name::kAligner,
                    enum_util::FormatEnumName(kBwa),
                    cli_opt_name::kAligner,
                    enum_util::FormatEnumName(kCustom),
                    cli_opt_name::kModel));
  }
  if (sample_type == kAutoDetect) {
    throw error::Error(
        fmt::format("Sample type '{}' is not supported for the tumor-only-wgs workflow. Use '{} {}' or '{} {}', or "
                    "'{} {}' with an explicit {}.",
                    enum_util::FormatEnumName(kAutoDetect),
                    cli_opt_name::kSampleType,
                    enum_util::FormatEnumName(kCellLine),
                    cli_opt_name::kSampleType,
                    enum_util::FormatEnumName(kFfpe),
                    cli_opt_name::kAligner,
                    enum_util::FormatEnumName(AlignerType::kCustom),
                    cli_opt_name::kModel));
  }
  const auto sample_type_str = enum_util::FormatEnumName(sample_type);
  const auto path =
      fmt::format("/resources/model-tumor-only-wgs-{}-{}.txt.gz", sample_type_str, AlignerToModelString(aligner));

  Logging::Info("Resolved model: {} based on aligner {} and sample_type {}",
                path,
                enum_util::FormatEnumName(aligner),
                sample_type_str);
  return fs::path{path};
}

namespace {
/**
 * @brief Resolve the paired germline SNV and indel model paths from filename templates.
 * @details Each template must contain two replacement fields, ordered as library type then aligner. Keeping this shared
 * implementation guarantees that the paired models use identical selection tokens.
 * @param library_type Concrete library type used in the bundled model filename.
 * @param aligner Aligner used in the bundled model filename.
 * @param snv_template Format template for the SNV model path.
 * @param indel_template Format template for the indel model path.
 * @return Resolved SNV and indel model paths, in that order.
 * @throws xoos::error::Error if the library type is unresolved or the aligner is custom.
 */
std::pair<fs::path, fs::path> ResolveGermlineModelPair(const LibraryType library_type,
                                                       const AlignerType aligner,
                                                       const std::string_view snv_template,
                                                       const std::string_view indel_template) {
  if (aligner == kCustom) {
    throw error::Error(fmt::format("Cannot resolve germline models for {} aligner — provide {} and {} explicitly",
                                   enum_util::FormatEnumName(aligner),
                                   cli_opt_name::kSnvModel,
                                   cli_opt_name::kIndelModel));
  }

  const auto library_type_str = LibraryTypeModelToken(library_type);
  const auto aligner_str = AlignerToModelString(aligner);

  const auto snv_path = fmt::format(fmt::runtime(snv_template), library_type_str, aligner_str);
  const auto indel_path = fmt::format(fmt::runtime(indel_template), library_type_str, aligner_str);

  Logging::Info("Resolved SNV model: {} for aligner {} and library type {}", snv_path, aligner_str, library_type_str);
  Logging::Info(
      "Resolved indel model: {} for aligner {} and library type {}", indel_path, aligner_str, library_type_str);

  return {fs::path{snv_path}, fs::path{indel_path}};
}
}  // namespace

std::pair<fs::path, fs::path> ResolveGermlineWgsModels(const LibraryType library_type,
                                                       const AlignerType aligner,
                                                       const GermlineModelType model_type) {
  using enum GermlineModelType;
  switch (model_type) {
    case kLarge:
      return ResolveGermlineModelPair(library_type,
                                      aligner,
                                      "/resources/model-germline-{}-{}-large-snv.txt.gz",
                                      "/resources/model-germline-{}-{}-large-indel.txt.gz");
    case kSmall:
      return ResolveGermlineModelPair(library_type,
                                      aligner,
                                      "/resources/model-germline-{}-{}-small-snv.txt.gz",
                                      "/resources/model-germline-{}-{}-small-indel.txt.gz");
    default:
      throw error::Error("Unknown germline WGS model type");
  }
}

std::string ModelThresholdsKey(const SampleType sample_type, const AlignerType aligner) {
  return fmt::format("{}-{}", enum_util::FormatEnumName(sample_type), enum_util::FormatEnumName(aligner));
}

std::optional<ModelThresholds> ResolveModelThresholds(const SampleType sample_type,
                                                      const AlignerType aligner,
                                                      const StrMap<ModelThresholds>& model_thresholds) {
  if (aligner == kCustom) {
    throw error::Error(fmt::format("Cannot resolve thresholds for {} aligner — provide explicit score thresholds",
                                   enum_util::FormatEnumName(aligner)));
  }

  if (model_thresholds.empty()) {
    return std::nullopt;
  }

  // Giraffe is rejected by the model resolvers before this is reached, so only bwa keys are resolved here.
  const auto key = ModelThresholdsKey(sample_type, aligner);
  const auto it = model_thresholds.find(key);
  if (it == model_thresholds.end()) {
    Logging::Warn("No model_thresholds entry for '{}', using top-level config defaults", key);
    return std::nullopt;
  }

  Logging::Info(
      "Resolved model thresholds for '{}': snv_min_ml_score={}, indel_min_ml_score={}, min_af={}, min_tumor_af={}",
      key,
      FormatOptionalThreshold(it->second.snv_min_ml_score),
      FormatOptionalThreshold(it->second.indel_min_ml_score),
      FormatOptionalThreshold(it->second.min_af),
      FormatOptionalThreshold(it->second.min_tumor_af));
  return it->second;
}

}  // namespace xoos::svc
