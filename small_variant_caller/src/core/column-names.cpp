#include "column-names.h"

namespace xoos::svc {

vec<std::string> FindUnsupportedFeatureNames(const vec<std::string>& names,
                                             const StrUnorderedSet& supported_names,
                                             const bool allow_tumor_normal_prefix) {
  vec<std::string> unsupported_names{};
  for (const auto& name : names) {
    if (supported_names.contains(name)) {
      continue;
    }
    if (allow_tumor_normal_prefix) {
      if (name.starts_with(kTumorPrefix) && supported_names.contains(name.substr(kTumorPrefix.size()))) {
        // feature for the tumor sample
        continue;
      }
      if (name.starts_with(kNormalPrefix) && supported_names.contains(name.substr(kNormalPrefix.size()))) {
        // feature for the normal sample
        continue;
      }
    }
    unsupported_names.emplace_back(name);
  }
  return unsupported_names;
}

vec<std::string> FindUnsupportedBamFeatureNames(const vec<std::string>& names, const bool allow_tumor_normal_prefix) {
  StrUnorderedSet supported_names = {};
  supported_names.insert(kVariantIDFeatureNames.begin(), kVariantIDFeatureNames.end());
  supported_names.insert(kBamFeatureNames.begin(), kBamFeatureNames.end());
  return FindUnsupportedFeatureNames(names, supported_names, allow_tumor_normal_prefix);
}

vec<std::string> FindUnsupportedVcfFeatureNames(const vec<std::string>& names) {
  StrUnorderedSet supported_names = {};
  supported_names.insert(kVariantIDFeatureNames.begin(), kVariantIDFeatureNames.end());
  supported_names.insert(kVcfFeatureNames.begin(), kVcfFeatureNames.end());
  return FindUnsupportedFeatureNames(names, supported_names, false);
}

vec<std::string> FindUnsupportedScoringFeatureNames(const vec<std::string>& names,
                                                    const bool allow_tumor_normal_prefix) {
  StrUnorderedSet supported_names = {};
  supported_names.insert(kVariantIDFeatureNames.begin(), kVariantIDFeatureNames.end());
  supported_names.insert(kBamFeatureNames.begin(), kBamFeatureNames.end());
  supported_names.insert(kVcfFeatureNames.begin(), kVcfFeatureNames.end());
  return FindUnsupportedFeatureNames(names, supported_names, allow_tumor_normal_prefix);
}

}  // namespace xoos::svc
