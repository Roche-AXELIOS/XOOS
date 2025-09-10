#include "column-names.h"

namespace xoos::svc {

/**
 * @brief Find feature names not supported.
 * @param names Vector of feature names to be evaluated
 * @param supported_names Set of feature names supported
 * @param allow_num_suffix Allow suffix of 1 or 2 in feature names
 * @return Vector of feature names not supported
 */
vec<std::string> FindUnsupportedFeatureNames(const vec<std::string>& names,
                                             const std::set<std::string, std::less<>>& supported_names,
                                             bool allow_num_suffix) {
  vec<std::string> unsupported_names{};
  for (const auto& name : names) {
    if (!supported_names.contains(name)) {
      if (allow_num_suffix) {
        if (name.ends_with("1") || name.ends_with("2")) {
          const std::string name_no_suffix = name.substr(0, name.size() - 1);
          if (!supported_names.contains(name_no_suffix)) {
            unsupported_names.emplace_back(name);
          }
        } else {
          unsupported_names.emplace_back(name);
        }
      } else {
        unsupported_names.emplace_back(name);
      }
    }
  }
  return unsupported_names;
}

/**
 * @brief Find BAM feature names not supported.
 * @param names Vector of feature names to be evaluated
 * @param allow_num_suffix Allow suffix of 1 or 2 in feature names
 * @return Vector of feature names not supported
 */
vec<std::string> FindUnsupportedBamFeatureNames(const vec<std::string>& names, bool allow_num_suffix) {
  StrSet supported_names = {};
  supported_names.insert(kVariantIDFeatureNames.begin(), kVariantIDFeatureNames.end());
  supported_names.insert(kBamFeatureNames.begin(), kBamFeatureNames.end());
  return FindUnsupportedFeatureNames(names, supported_names, allow_num_suffix);
}

/**
 * @brief Find VCF feature names not supported.
 * @param names Vector of feature names to be evaluated
 * @param allow_num_suffix Allow suffix of 1 or 2 in feature names
 * @return Vector of feature names not supported
 */
vec<std::string> FindUnsupportedVcfFeatureNames(const vec<std::string>& names, bool allow_num_suffix) {
  StrSet supported_names = {};
  supported_names.insert(kVariantIDFeatureNames.begin(), kVariantIDFeatureNames.end());
  supported_names.insert(kVcfFeatureNames.begin(), kVcfFeatureNames.end());
  return FindUnsupportedFeatureNames(names, supported_names, allow_num_suffix);
}

/**
 * @brief Find scoring feature names not supported.
 * @param names Vector of feature names to be evaluated
 * @param allow_num_suffix Allow suffix of 1 or 2 in feature names
 * @return Vector of feature names not supported
 */
vec<std::string> FindUnsupportedScoringFeatureNames(const vec<std::string>& names, bool allow_num_suffix) {
  StrSet supported_names = {};
  supported_names.insert(kVariantIDFeatureNames.begin(), kVariantIDFeatureNames.end());
  supported_names.insert(kBamFeatureNames.begin(), kBamFeatureNames.end());
  supported_names.insert(kVcfFeatureNames.begin(), kVcfFeatureNames.end());
  return FindUnsupportedFeatureNames(names, supported_names, allow_num_suffix);
}

}  // namespace xoos::svc
