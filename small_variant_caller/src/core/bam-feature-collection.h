#pragma once

#include <optional>

#include "variant-info.h"

namespace xoos::svc {

/**
 * @brief Tuple of BAM features including variant and reference allele features.
 */
struct BamFeatureTuple {
  // Variant allele features
  UnifiedVariantFeature var_feat{};
  // Reference allele features
  UnifiedReferenceFeature ref_feat{};
};

/**
 * @brief Tuple of BAM features including tumor-specific and normal-specific variant and reference allele features.
 */
struct TumorNormalBamFeatureTuple : BamFeatureTuple {
  // Tumor-specific variant allele features
  UnifiedVariantFeature tumor_var_feat{};
  // Tumor-specific reference allele features
  UnifiedReferenceFeature tumor_ref_feat{};
  // Normal-specific variant allele features
  UnifiedVariantFeature normal_var_feat{};
  // Normal-specific reference allele features
  UnifiedReferenceFeature normal_ref_feat{};
};

/**
 * @brief Collection of extracted BAM features for a given chromosomal region including variant and reference allele
 * features for overall, tumor-specific, and normal-specific contexts.
 * @note All features are interpreted as being from the same chromosomal region. There, no chromosome information is
 * stored here, and the reference allele features are indexed only by position.
 */
struct BamRegionFeatureCollection {
 public:
  // Variant features indexed by Variant ID
  VarIdToVarBamFeatures var_features{};
  // Reference allele features indexed by chromosomal position
  PosToRefBamFeatures ref_features{};
  // Tumor-specific variant features indexed by Variant ID
  VarIdToVarBamFeatures tumor_var_features{};
  // Tumor-specific reference allele features indexed by chromosomal position
  PosToRefBamFeatures tumor_ref_features{};
  // Normal-specific variant features indexed by Variant ID
  VarIdToVarBamFeatures normal_var_features{};
  // Normal-specific reference allele features indexed by chromosomal position
  PosToRefBamFeatures normal_ref_features{};

  /**
   * @brief Retrieves the BAM feature tuple for the specified Variant ID.
   * @param vid Variant ID for which to retrieve features
   * @return BAM feature tuple if variant is found, otherwise std::nullopt
   */
  std::optional<BamFeatureTuple> GetBamFeatureTuple(const VariantId& vid) const;

  /**
   * @brief Retrieves the tumor-normal BAM feature tuple for the specified Variant ID.
   * @param vid Variant ID for which to retrieve features
   * @return Tumor-normal BAM feature tuple if variant found, otherwise std::nullopt
   */
  std::optional<TumorNormalBamFeatureTuple> GetTumorNormalBamFeatureTuple(const VariantId& vid) const;

  /**
   * @brief Clears all stored BAM features from the collection.
   */
  void Clear();
};

/**
 * @brief Find the variant allele BAM feature for a variant.
 * @param features Map of variant allele BAM features indexed by Variant ID
 * @param vid Variant ID
 * @return UnifiedVariantFeature for the VariantId or a zeroed feature if not found
 */
const UnifiedVariantFeature& FindVarBamFeature(const VarIdToVarBamFeatures& features, const VariantId& vid);

/**
 * @brief Find the reference allele BAM feature for a variant.
 * @param features Map of reference allele features indexed by chromosomal position
 * @param vid Variant ID
 * @return UnifiedReferenceFeature for the VariantId or a zeroed feature if not found
 */
const UnifiedReferenceFeature& FindRefBamFeature(const PosToRefBamFeatures& features, const VariantId& vid);

/**
 * @brief Find the reference allele BAM feature for a variant.
 * @param features Nested map of reference allele features indexed by chromosome name and position
 * @param vid Variant ID
 * @return UnifiedReferenceFeature for the VariantId or a zeroed feature if not found
 */
const UnifiedReferenceFeature& FindRefBamFeature(const ChromPosToRefBamFeatures& features, const VariantId& vid);

}  // namespace xoos::svc
