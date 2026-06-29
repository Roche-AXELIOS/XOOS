#pragma once

#include <xoos/types/vec.h>

#include "bam-feature-collection.h"
#include "feature-normalization.h"
#include "variant-info.h"
#include "xoos/types/float.h"

namespace xoos::svc {

/**
 * @brief Generate a vector of feature values for a given variant.
 * @param feature_cols Vector of feature columns
 * @param vid Variant ID
 * @param bam_feat Variant's BAM features in general
 * @param vcf_feat Variant's VCF features
 * @param normalize_target Target depth for feature value normalization
 * @return Vector of feature values
 */
vec<f64> GetFeatureVec(const vec<FeatureColumn>& feature_cols,
                       const VariantId& vid,
                       const BamFeatureTuple& bam_feat,
                       const VcfFeature& vcf_feat,
                       const DepthTuple& normalize_target);

/**
 * @brief Generate a vector of feature values for a given variant.
 * @param feature_cols Vector of feature columns
 * @param vid Variant ID
 * @param bam_feat Variant's BAM features in tumor-normal sample context
 * @param vcf_feat Variant's VCF features
 * @param normalize_target Target depth for feature value normalization
 * @return Vector of feature values
 */
vec<f64> GetFeatureVec(const vec<FeatureColumn>& feature_cols,
                       const VariantId& vid,
                       const TumorNormalBamFeatureTuple& bam_feat,
                       const VcfFeature& vcf_feat,
                       const DepthTuple& normalize_target);

}  // namespace xoos::svc
