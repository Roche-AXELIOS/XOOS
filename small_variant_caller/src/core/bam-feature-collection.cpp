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

namespace {

// Look up the reference-feature bucket holding per-position aggregate features. These are keyed by the raw variant
// position (GetAggFeaturePos), which differs from GetRefFeaturePos for deletions.
const UnifiedReferenceFeature& FindAggRefBamFeature(const PosToRefBamFeatures& features, const VariantId& vid) {
  const auto it = features.find(vid.GetAggFeaturePos());
  if (it == features.end()) {
    return kZeroUnifiedReferenceFeature;
  }
  return it->second;
}

// Overlay the per-position aggregate feature fields from the raw-position bucket onto the reference feature returned
// for a variant. For SNVs and insertions the two buckets coincide (GetRefFeaturePos == GetAggFeaturePos), so this is a
// no-op copy; for deletions it pulls the aggregate values from the raw variant position.
void OverlayAggregateFeatures(UnifiedReferenceFeature& dst, const UnifiedReferenceFeature& agg) {
  dst.del_1bp_support = agg.del_1bp_support;
  dst.del_2bp_support = agg.del_2bp_support;
  dst.del_3bp_support = agg.del_3bp_support;
  dst.del_gt3bp_support = agg.del_gt3bp_support;
  dst.ins_1bp_support = agg.ins_1bp_support;
  dst.ins_2bp_support = agg.ins_2bp_support;
  dst.ins_3bp_support = agg.ins_3bp_support;
  dst.ins_gt3bp_support = agg.ins_gt3bp_support;
  dst.del_1bp_af = agg.del_1bp_af;
  dst.del_2bp_af = agg.del_2bp_af;
  dst.del_3bp_af = agg.del_3bp_af;
  dst.del_gt3bp_af = agg.del_gt3bp_af;
  dst.ins_1bp_af = agg.ins_1bp_af;
  dst.ins_2bp_af = agg.ins_2bp_af;
  dst.ins_3bp_af = agg.ins_3bp_af;
  dst.ins_gt3bp_af = agg.ins_gt3bp_af;
  dst.del_1bp_ins_1bp_af_logratio = agg.del_1bp_ins_1bp_af_logratio;
  dst.del_1bp_2bp_af_logratio = agg.del_1bp_2bp_af_logratio;
  dst.del_2bp_3bp_af_logratio = agg.del_2bp_3bp_af_logratio;
  dst.del_3bp_gt3bp_af_logratio = agg.del_3bp_gt3bp_af_logratio;
  dst.ins_1bp_2bp_af_logratio = agg.ins_1bp_2bp_af_logratio;
  dst.ins_2bp_3bp_af_logratio = agg.ins_2bp_3bp_af_logratio;
  dst.ins_3bp_gt3bp_af_logratio = agg.ins_3bp_gt3bp_af_logratio;
  dst.snv_support = agg.snv_support;
  dst.snv_af = agg.snv_af;
}

}  // namespace

std::optional<BamFeatureTuple> BamRegionFeatureCollection::GetBamFeatureTuple(const VariantId& vid) const {
  const auto it = var_features.find(vid);
  if (it == var_features.end()) {
    return std::nullopt;
  }
  // Variant found, construct the feature tuple
  BamFeatureTuple tuple{};
  tuple.var_feat = it->second;
  tuple.ref_feat = FindRefBamFeature(ref_features, vid);
  // Per-position aggregate features are keyed by the raw variant position, which differs from the reference-feature
  // position for deletions. Overlay them from the raw-position bucket.
  OverlayAggregateFeatures(tuple.ref_feat, FindAggRefBamFeature(ref_features, vid));
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
  // Per-position aggregate features are keyed by the raw variant position, which differs from the reference-feature
  // position for deletions. Overlay them from the raw-position bucket for each sample context.
  OverlayAggregateFeatures(tuple.ref_feat, FindAggRefBamFeature(ref_features, vid));
  OverlayAggregateFeatures(tuple.tumor_ref_feat, FindAggRefBamFeature(tumor_ref_features, vid));
  OverlayAggregateFeatures(tuple.normal_ref_feat, FindAggRefBamFeature(normal_ref_features, vid));
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
