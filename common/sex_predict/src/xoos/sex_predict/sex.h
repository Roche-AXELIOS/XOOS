#pragma once

#include <string>

namespace xoos::sex_predict {

enum class Sex {
  kUnknown,
  kMale,
  kFemale
};
enum class SexChromHandling {
  kAllChrom,
  kNoSexChrom,
  kOnlySexChrom,
  kNoChromY
};

Sex ParseSexOption(const std::string& sex);
std::string GetDescriptionForSex(Sex sex);

}  // namespace xoos::sex_predict
