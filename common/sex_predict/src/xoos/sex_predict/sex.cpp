#include "sex.h"

#include <algorithm>
#include <array>
#include <string>
#include <string_view>

#include <xoos/genome/region.h>

namespace xoos::sex_predict {

std::string GetDescriptionForSex(const Sex sex) {
  switch (sex) {
    case Sex::kMale:
      return "Male";
    case Sex::kFemale:
      return "Female";
    default:
      return "Unknown";
  }
}

bool IsInAllosome(const std::string_view region) {
  constexpr std::array<std::string_view, 4> kAllosomes = {"chrX", "chrY", "X", "Y"};
  return std::ranges::any_of(
      kAllosomes, [region](const std::string_view allosome) { return genome::IsEqualContig(region, allosome); });
}

bool IsInChromY(const std::string_view region) {
  constexpr std::array<std::string_view, 2> kChromYNames = {"chrY", "Y"};
  return std::ranges::any_of(kChromYNames,
                             [region](const std::string_view name) { return genome::IsEqualContig(region, name); });
}

}  // namespace xoos::sex_predict
