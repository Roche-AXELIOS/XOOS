#include "bam-feature-collection.h"

namespace xoos::svc {

const UnifiedVariantFeature& FindVarBamFeature(const VarIdToVarBamFeatures& features, const VariantId& vid) {
  const auto it = features.find(vid);
  if (it == features.end()) {
    return kZeroUnifiedVariantFeature;
  }
  return it->second;
}

const UnifiedReferenceFeature& FindRefBamFeature(const PosToRefBamFeatures& features, const VariantId& vid) {
  const auto it = features.find(vid.GetRefFeaturePos());
  if (it == features.end()) {
    return kZeroUnifiedReferenceFeature;
  }
  return it->second;
}

const UnifiedReferenceFeature& FindRefBamFeature(const ChromPosToRefBamFeatures& features, const VariantId& vid) {
  const auto it = features.find(vid.chrom);
  if (it == features.end()) {
    return kZeroUnifiedReferenceFeature;
  }
  return FindRefBamFeature(it->second, vid);
}

std::optional<BamFeatureTuple> BamRegionFeatureCollection::GetBamFeatureTuple(const VariantId& vid) const {
  const auto it = var_features.find(vid);
  if (it == var_features.end()) {
    return std::nullopt;
  }
  // Variant found, construct the feature tuple
  BamFeatureTuple tuple{};
  tuple.var_feat = it->second;
  tuple.ref_feat = FindRefBamFeature(ref_features, vid);
  return tuple;
}

std::optional<TumorNormalBamFeatureTuple> BamRegionFeatureCollection::GetTumorNormalBamFeatureTuple(
    const VariantId& vid) const {
  const auto it = var_features.find(vid);
  if (it == var_features.end()) {
    return std::nullopt;
  }
  // Variant found, construct the feature tuple
  TumorNormalBamFeatureTuple tuple{};
  tuple.var_feat = it->second;
  tuple.ref_feat = FindRefBamFeature(ref_features, vid);
  tuple.tumor_var_feat = FindVarBamFeature(tumor_var_features, vid);
  tuple.tumor_ref_feat = FindRefBamFeature(tumor_ref_features, vid);
  tuple.normal_var_feat = FindVarBamFeature(normal_var_features, vid);
  tuple.normal_ref_feat = FindRefBamFeature(normal_ref_features, vid);
  return tuple;
}

void BamRegionFeatureCollection::Clear() {
  var_features.clear();
  ref_features.clear();
  tumor_var_features.clear();
  tumor_ref_features.clear();
  normal_var_features.clear();
  normal_ref_features.clear();
}

}  // namespace xoos::svc
