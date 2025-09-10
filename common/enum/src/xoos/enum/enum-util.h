#pragma once

#include <algorithm>
#include <optional>
#include <string>
#include <vector>

#include <fmt/ranges.h>

#include <magic_enum/magic_enum.hpp>

#include <xoos/log/logging.h>

namespace xoos::enum_util {

/**
 * @brief Parses a string into an enum value. The string can be the enum name or the enum name without the leading 'k'.
 */
template <typename T>
std::optional<T> ParseEnumName(const std::string& text, const std::optional<T>& default_value = std::nullopt) {
  auto k_text = "k" + text;
  std::erase_if(k_text, [](char c) { return c == '-' || c == '.'; });
  auto result = magic_enum::enum_cast<T>(k_text, magic_enum::case_insensitive);
  return result ? result : default_value;
}

/**
 * @brief Inserts a dash before each uppercase letter except for the first letter and makes first letter lowercase
 */
std::string FormatEnumNameLowerDash(const std::string& name);

/**
 * @brief Converts enum values into string name
 */
template <typename T>
std::string FormatEnumName(T value) {
  auto name = std::string{magic_enum::enum_name(value)};
  name = name.starts_with("k") ? std::string{name.substr(1)} : name;
  name = FormatEnumNameLowerDash(name);
  return name;
}

/**
 * @brief Converts vector of enum values into string name
 */
template <typename T>
std::vector<std::string> FormatEnumName(const std::vector<T>& values) {
  std::vector<std::string> names;
  names.reserve(values.size());
  for (const auto& value : values) {
    names.emplace_back(FormatEnumName(value));
  }
  return names;
}

/**
 * @brief Formats all enum names into a string.
 */
template <typename T>
std::string FormatEnumNames(const std::string& l_wrap = "[", const std::string& r_wrap = "]") {
  auto values = magic_enum::enum_values<T>();
  std::vector<std::string> names{values.size()};
  std::ranges::transform(
      std::cbegin(values), std::cend(values), std::begin(names), [](const auto& v) { return FormatEnumName(v); });
  return fmt::format("{}{}{}", l_wrap, fmt::join(names, ", "), r_wrap);
}

/**
 * @brief Parses cli into an enum value.
 */
template <typename T>
std::optional<T> ParseEnumNameCliOpt(const std::string& cli_opt_name, const std::string& text) {
  auto value = ParseEnumName<T>(text);
  if (!value) {
    Logging::Warn("Invalid value for '{}'. Must be one of {}", cli_opt_name, FormatEnumNames<T>());
  }
  return value;
}

/**
 * @brief Parses a vector of strings into a vector of enum values.
 */
template <typename T>
std::optional<std::vector<T>> ParseEnumNameCliOpt(const std::string& cli_opt_name,
                                                  const std::vector<std::string>& texts) {
  std::vector<T> values;
  values.reserve(texts.size());
  for (const auto& text : texts) {
    auto value = ParseEnumNameCliOpt<T>(cli_opt_name, text);
    if (!value) {
      return std::nullopt;
    }
    values.emplace_back(*value);
  }
  return values;
}

}  // namespace xoos::enum_util
