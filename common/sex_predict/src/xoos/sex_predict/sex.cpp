#include "sex.h"

#include <map>
#include <vector>

#include <xoos/log/logging.h>
#include <xoos/util/string-functions.h>

namespace xoos::sex_predict {
const std::vector kPossibleSexChars{'M', 'F', 'N'};
const std::map<char, Sex> kCharToSex{{'M', Sex::kMale}, {'F', Sex::kFemale}, {'N', Sex::kUnknown}};

Sex ParseSexOption(const std::string& sex) {
  if (sex.empty()) {
    Logging::Info("No sex specified, defaulting to 'N'.");
    return Sex::kUnknown;
  }
  auto choose_sex = [&](char sex_char) -> Sex {
    const auto it = kCharToSex.find(sex_char);
    if (it == kCharToSex.end()) {
      throw std::runtime_error(fmt::format("invalid sex character: {}", sex_char));
    }
    return it->second;
  };
  // set the first character to upper_case
  // extract the first character of sex and set it to upper case
  return choose_sex(static_cast<char>(std::toupper(sex[0])));
}

std::string GetDescriptionForSex(const Sex sex) {
  // get value for sex from kCharToSex
  switch (sex) {
    case Sex::kMale:
      return "Male";
    case Sex::kFemale:
      return "Female";
    default:
      return "Unknown";
  }
}
}  // namespace xoos::sex_predict
