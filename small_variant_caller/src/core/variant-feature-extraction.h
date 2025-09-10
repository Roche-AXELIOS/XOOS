#pragma once

#include <optional>

#include <xoos/types/vec.h>

#include "variant-info.h"

namespace xoos::svc {

vec<double> GetFeatureVec(const vec<UnifiedFeatureCols>& feature_names,
                          const VariantId& vid,
                          const UnifiedVariantFeature& bam_feat,
                          const UnifiedReferenceFeature& ref_feat,
                          const VcfFeature& vcf_feat,
                          std::optional<u32> normalize_target);
size_t GetFeatureVecLength(const vec<UnifiedFeatureCols>& feature_names);

}  // namespace xoos::svc
