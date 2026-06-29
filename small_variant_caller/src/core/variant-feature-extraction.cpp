#include "variant-feature-extraction.h"

#include <cmath>

#include <xoos/error/error.h>

#include "variant-info-serializer.h"

namespace xoos::svc {

vec<f64> GetFeatureVec(const vec<FeatureColumn>& feature_cols,
                       const VariantId& vid,
                       const BamFeatureTuple& bam_feat,
                       const VcfFeature& vcf_feat,
                       const DepthTuple& normalize_target) {
  using enum SampleContext;
  vec<f64> feature_vals;
  feature_vals.reserve(feature_cols.size());
  for (const auto& col : feature_cols) {
    f64 val =
        VariantInfoSerializer::NumericalizeFeature(col.enum_val, vid, bam_feat.var_feat, vcf_feat, bam_feat.ref_feat);
    if (std::isnan(val)) {
      const auto& serialized = VariantInfoSerializer::SerializeBamFeatureRow({col}, vid, bam_feat);
      throw error::Error("Requested Field '{}' Value '{}' cannot be numericalized for Variant '{}'",
                         GetFeatureName(col),
                         serialized.front(),
                         vid.ToString());
    }
    if (normalize_target.total > 0 && kNormalizableFeatureCols.contains(col.enum_val)) {
      feature_vals.push_back(val / normalize_target.total);
    } else {
      feature_vals.push_back(val);
    }
  }
  return feature_vals;
}

vec<f64> GetFeatureVec(const vec<FeatureColumn>& feature_cols,
                       const VariantId& vid,
                       const TumorNormalBamFeatureTuple& bam_feat,
                       const VcfFeature& vcf_feat,
                       const DepthTuple& normalize_target) {
  using enum SampleContext;
  vec<f64> feature_vals;
  feature_vals.reserve(feature_cols.size());
  for (const auto& col : feature_cols) {
    f64 val = 0;
    switch (col.sample_context) {
      case kTumor:
        val = VariantInfoSerializer::NumericalizeFeature(
            col.enum_val, vid, bam_feat.tumor_var_feat, vcf_feat, bam_feat.tumor_ref_feat);
        break;
      case kNormal:
        val = VariantInfoSerializer::NumericalizeFeature(
            col.enum_val, vid, bam_feat.normal_var_feat, vcf_feat, bam_feat.normal_ref_feat);
        break;
      case kNone:
      default:
        val = VariantInfoSerializer::NumericalizeFeature(
            col.enum_val, vid, bam_feat.var_feat, vcf_feat, bam_feat.ref_feat);
        break;
    }
    if (std::isnan(val)) {
      const auto& serialized = VariantInfoSerializer::SerializeBamFeatureRow({col}, vid, bam_feat);
      throw error::Error("Requested Field '{}' Value '{}' cannot be numericalized for Variant '{}'",
                         GetFeatureName(col),
                         serialized.front(),
                         vid.ToString());
    }
    if (kNormalizableFeatureCols.contains(col.enum_val)) {
      // normalize the feature value based on the feature's sample context, if applicable
      // feature value will not be normalized if the normalization target for the feature's sample context is zero
      if (col.sample_context == kTumor && normalize_target.tumor > 0) {
        feature_vals.push_back(val / normalize_target.tumor);
      } else if (col.sample_context == kNormal && normalize_target.normal > 0) {
        feature_vals.push_back(val / normalize_target.normal);
      } else if (col.sample_context == kNone && normalize_target.total > 0) {
        feature_vals.push_back(val / normalize_target.total);
      } else {
        feature_vals.push_back(val);
      }
    } else {
      feature_vals.push_back(val);
    }
  }
  return feature_vals;
}

}  // namespace xoos::svc
