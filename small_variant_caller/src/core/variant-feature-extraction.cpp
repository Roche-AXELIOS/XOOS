#include "variant-feature-extraction.h"

#include <cmath>

#include <xoos/error/error.h>

#include "core/variant-info-serializer.h"

namespace xoos::svc {

// dummy structs
static const VariantId kDummyVid = VariantId();
static const UnifiedVariantFeature kDummyBamFeat = UnifiedVariantFeature();
static const VcfFeature kDummyVcfFeat = VcfFeature();
static const UnifiedReferenceFeature kDummyRefFeat = UnifiedReferenceFeature();

/**
 * @brief Generate a vector of feature values for a given variant.
 * @param feature_names Vector of feature names
 * @param vid Variant ID
 * @param bam_feat Variant's BAM features
 * @param ref_feat Variant's reference features
 * @param vcf_feat Variant's VCF features
 * @param normalize_target Target depth for feature value normalization
 * @return Vector of feature values
 */
vec<double> GetFeatureVec(const vec<UnifiedFeatureCols>& feature_names,
                          const VariantId& vid,
                          const UnifiedVariantFeature& bam_feat,
                          const UnifiedReferenceFeature& ref_feat,
                          const VcfFeature& vcf_feat,
                          std::optional<u32> normalize_target) {
  // cannot normalize feature values to target depth of 0
  const bool normalize{normalize_target.has_value() && normalize_target.value() > 0};
  vec<double> feature_vals;
  feature_vals.reserve(feature_names.size());
  for (const auto& col : feature_names) {
    auto val = VariantInfoSerializer::NumericalizeFeature(col, vid, bam_feat, vcf_feat, ref_feat);
    if (std::isnan(val)) {
      throw error::Error("Requested Field " + GetFeatureName(col) + " cannot be numericalized");
    }
    if (normalize && kNormalizableFeatureCols.contains(col)) {
      feature_vals.push_back(val / normalize_target.value());
    } else {
      feature_vals.push_back(val);
    }
  }
  return feature_vals;
}

/**
 * @brief Return the length of feature vector.
 * @param feature_names Feature names
 * @return Length of feature vector
 */
size_t GetFeatureVecLength(const vec<UnifiedFeatureCols>& feature_names) {
  return GetFeatureVec(feature_names, kDummyVid, kDummyBamFeat, kDummyRefFeat, kDummyVcfFeat, 0).size();
}

}  // namespace xoos::svc
