#pragma once

#include <optional>
#include <utility>

#include <xoos/types/fs.h>

#include "core/aligner-type.h"
#include "core/config.h"
#include "core/germline-model-type.h"
#include "core/library-type.h"
#include "core/sample-type.h"

namespace xoos::svc {

/**
 * @brief Resolve the tumor-normal-wgs model path from sample type and aligner type.
 * @details Maps the combination of sample type (ffpe, cell-line) and aligner (bwa) to the corresponding
 *          model file in /resources/. Giraffe is not supported for tumor-normal-wgs because no pre-trained
 *          giraffe model exists; pass kCustom with an explicit --model for giraffe-aligned data.
 * @param sample_type The sample preparation type.
 * @param aligner The aligner type. Must not be kGiraffe (unsupported for tumor-normal-wgs) or kCustom
 *                (custom requires explicit --model).
 * @return Path to the model file (e.g., "/resources/model-tumor-normal-wgs-ffpe-bwa.txt.gz").
 * @throws xoos::Error if aligner is kGiraffe or kCustom.
 */
fs::path ResolveTumorNormalModel(SampleType sample_type, AlignerType aligner);

/**
 * @brief Resolve the tumor-only-wgs model path from sample type and aligner type.
 * @details Maps the combination of sample type (ffpe, cell-line) and aligner (bwa) to the corresponding
 *          model file in /resources/. Giraffe is not supported because no pre-trained giraffe model exists;
 *          pass kCustom with an explicit --model for giraffe-aligned data.
 * @param sample_type The sample preparation type.
 * @param aligner The aligner type. Must not be kGiraffe (unsupported for tumor-only-wgs) or kCustom
 *                (custom requires explicit --model).
 * @return Path to the model file (e.g., "/resources/model-tumor-only-wgs-cell-line-bwa.txt.gz").
 * @throws xoos::Error if aligner is kGiraffe or kCustom, or if sample_type is kAutoDetect.
 */
fs::path ResolveTumorOnlyWgsModel(SampleType sample_type, AlignerType aligner);

/**
 * @brief Resolve germline WGS SNV and indel model paths from library type, aligner, and model family.
 * @details Maps the concrete library type, aligner, and germline model type to the corresponding bundled model pair.
 * @param library_type Concrete library type. Auto-detect must be resolved before calling.
 * @param aligner The aligner type. Must not be kCustom (custom requires explicit --snv-model/--indel-model).
 * @param model_type Bundled model family to select.
 * @return Pair of (SNV model path, indel model path).
 * @throws xoos::Error if the selection has no bundled model.
 */
std::pair<fs::path, fs::path> ResolveGermlineWgsModels(LibraryType library_type,
                                                       AlignerType aligner,
                                                       GermlineModelType model_type);

/**
 * @brief Build a model-threshold lookup key from sample type and aligner.
 * @details Produces the workflow-agnostic key used by tumor-normal-wgs and tumor-only-wgs model threshold
 * configuration. ResolveModelThresholds() uses the same key at runtime.
 * @param sample_type Sample preparation type, such as SampleType::kFfpe or SampleType::kCellLine.
 * @param aligner Aligner type, such as AlignerType::kBwa. Custom aligners are expected to use explicit thresholds.
 * @return String key in the format "{sample_type}-{aligner}".
 */
std::string ModelThresholdsKey(SampleType sample_type, AlignerType aligner);

/**
 * @brief Resolve per-model score thresholds from configuration.
 * @details Looks up the threshold entry keyed by ModelThresholdsKey(). Returns no value when the map is
 * empty or the key is absent, allowing callers to keep top-level defaults. Custom aligners require explicit
 * score thresholds and are rejected.
 * @param sample_type Sample preparation type.
 * @param aligner Aligner type. Must not be AlignerType::kCustom.
 * @param model_thresholds Per-model threshold map from SVCConfig.
 * @return Optional ModelThresholds for the sample-type/aligner combination.
 * @throws xoos::error::Error when @p aligner is AlignerType::kCustom.
 */
std::optional<ModelThresholds> ResolveModelThresholds(SampleType sample_type,
                                                      AlignerType aligner,
                                                      const StrMap<ModelThresholds>& model_thresholds);

}  // namespace xoos::svc
