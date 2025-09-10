#include "xoos/enum/enum-util.h"

#include <regex>

namespace xoos::enum_util {

/**
 * @brief Splits a string where the letter case changes, when it changes from upper to lower. Examples:
 * "HelloWorld" -> {"Hello", "World"}
 * "HelloWorld123" -> {"Hello", "World123"}
 * "Hello123World" -> {"Hello123", "World"}
 * "HelloWORLd" -> {"Hello", "WOR", "Ld"}
 */
static std::vector<std::string> SplitOnLetterChangeCase(const std::string& input) {
  const std::regex re("([A-Z]+(?![a-z])(?![\\d])|[A-Z][a-z]*[\\d]*)");
  const std::sregex_token_iterator first{input.begin(), input.end(), re};
  return {first, std::sregex_token_iterator{}};
}

/**
 * @brief Checks if the input is an acronym, an acronym is defined as a string with more than 1
 * letters and where all characters are uppercase
 */
static bool IsAcronym(const std::string_view& input) {
  return input.size() > 1 && std::ranges::all_of(input.begin(), input.end(), [](char c) { return std::isupper(c); });
}

/**
 * @brief changes the first letter to lowercase if it is not an acronym
 */
static std::string MakeWordLowerCaseIfNotAcronym(const std::string& input) {
  if (IsAcronym(input)) {
    return input;
  }
  std::string lower_input{input};
  lower_input[0] = static_cast<char>(std::tolower(input[0]));
  return lower_input;
}

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
