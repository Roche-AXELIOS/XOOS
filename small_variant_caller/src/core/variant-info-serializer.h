#pragma once

#include <map>
#include <string>
#include <utility>

#include "core/variant-info.h"

namespace xoos::svc {

/**
 * Synopsis:
 * This file provides serialization and deserialization methods for variant features,
 * reference features, and variant IDs. It also includes methods to load features from files
 * and numericalize feature values.
 */

class VariantInfoSerializer {
 public:
  static UnifiedVariantFeature DeserializeVariantFeature(const std::map<UnifiedFeatureCols, std::string>& serialized);
  static UnifiedReferenceFeature DeserializeRefFeature(const std::map<UnifiedFeatureCols, std::string>& serialized);
  static VariantId DeserializeVariantId(const std::map<UnifiedFeatureCols, std::string>& serialized);
  static std::pair<ChromToVariantInfoMap, RefInfoMap> LoadFeatures(const std::string& features_file);
  static std::pair<ChromToVariantInfoMap, RefInfoMap> LoadFeatures(const std::string& features_file,
                                                                   const ChromToVcfFeaturesMap& chr_to_vcf_feat_map);
  static ChromToVcfFeaturesMap LoadVcfFeatures(const std::string& features_file);
  static std::map<UnifiedFeatureCols, std::string> SerializeVariantFeature(const UnifiedVariantFeature& variant_info);
  static std::map<UnifiedFeatureCols, std::string> SerializeReferenceFeature(const UnifiedReferenceFeature& ref_info);
  static std::map<UnifiedFeatureCols, std::string> SerializeVariantId(const VariantId& vid);
  static double NumericalizeFeature(UnifiedFeatureCols col,
                                    const VariantId& vid,
                                    const UnifiedVariantFeature& bam_feat,
                                    const VcfFeature& vcf_feat,
                                    const UnifiedReferenceFeature& ref_feat);
  static double NumericalizeFeature(UnifiedFeatureCols col, const VariantId& vid);
  static double NumericalizeFeature(UnifiedFeatureCols col, const UnifiedVariantFeature& feat);
  static double NumericalizeFeature(UnifiedFeatureCols col, const VcfFeature& feat);
  static double NumericalizeFeature(UnifiedFeatureCols col, const UnifiedReferenceFeature& feat);
  static VcfFeature DeserializeVcfFeature(std::map<UnifiedFeatureCols, std::string>& serialized);
  static std::map<UnifiedFeatureCols, std::string> SerializeVcfFeature(const VcfFeature& feature);
};

}  // namespace xoos::svc
