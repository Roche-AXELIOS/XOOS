#include "xoos/enum/enum-util.h"

#include <algorithm>
#include <cctype>
#include <locale>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace xoos::enum_util {

namespace {

/**
 * @brief Splits a string where the letter case changes, when it changes from upper to lower. Examples:
 * "HelloWorld" -> {"Hello", "World"}
 * "HelloWorld123" -> {"Hello", "World123"}
 * "Hello123World" -> {"Hello123", "World"}
 * "HelloWORLd" -> {"Hello", "WOR", "Ld"}
 */
std::vector<std::string> SplitOnLetterChangeCase(const std::string_view input) {
  static const std::regex kPattern("([A-Z]+(?![a-z])(?![\\d])|[A-Z][a-z]*[\\d]*)");
  const std::cregex_token_iterator first{input.data(), input.data() + input.size(), kPattern};
  return {first, std::cregex_token_iterator{}};
}

/**
 * @brief Checks if the input is an acronym, an acronym is defined as a string with more than 1
 * letters and where all characters are uppercase
 */
bool IsAcronym(const std::string_view input) {
  return input.size() > 1 &&
         std::ranges::all_of(input, [](const char c) { return std::isupper(c, std::locale::classic()); });
}

/**
 * @brief changes the first letter to lowercase if it is not an acronym
 */
std::string MakeWordLowerCaseIfNotAcronym(const std::string_view input) {
  if (IsAcronym(input)) {
    return std::string{input};
  }
  std::string lower_input{input};
  lower_input[0] = std::tolower(input[0], std::locale::classic());
  return lower_input;
}

}  // namespace

/**
 * @brief Inserts a dash before each capital letter except for the first letter and makes first letter lowercase after
 * each dash if the word is not an acronym
 */
std::string FormatEnumNameLowerDash(const std::string& name) {
  const std::vector<std::string> parts = SplitOnLetterChangeCase(name);
  auto it = std::cbegin(parts);

  std::ostringstream formatted_name;
  if (it != std::cend(parts)) {
    formatted_name << MakeWordLowerCaseIfNotAcronym(*it);
    ++it;
  }
  for (; it != std::cend(parts); ++it) {
    formatted_name << "-" << MakeWordLowerCaseIfNotAcronym(*it);
  }
  return formatted_name.str();
}

}  // namespace xoos::enum_util
