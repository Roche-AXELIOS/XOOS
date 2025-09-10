#pragma once

#include <string>
#include <vector>

#include <fmt/format.h>

#include <CLI/CLI.hpp>

#include <xoos/enum/enum-util.h>

#include "xoos/cli/cli.h"

namespace xoos::cli {

using enum_util::FormatEnumName;
using enum_util::FormatEnumNames;
using enum_util::ParseEnumName;

template <class T>
T ParseEnumNameOrThrow(const std::string& name, const std::string& input) {
  auto input_parsed = ParseEnumName<T>(input);
  if (!input_parsed) {
    // Create a validation error message that has the same structure as CLI::CheckedTransformer
    auto msg = fmt::format("{}: Check {} value in {} FAILED", name, input, FormatEnumNames<T>("{", "}"));
    throw CLI::ValidationError(msg);
  }
  return *input_parsed;
}

template <class T>
std::function<void(const std::string&)> EnumOptionFunction(const std::string& name, T& value) {
  return [name, &value](const std::string& input) { value = ParseEnumNameOrThrow<T>(name, input); };
}

template <class T>
std::function<void(const std::vector<std::string>&)> EnumOptionFunction(const std::string& name,
                                                                        std::vector<T>& values) {
  return [name, &values](const std::vector<std::string>& inputs) {
    std::vector<T> inputs_parsed;
    for (const auto& input : inputs) {
      inputs_parsed.push_back(ParseEnumNameOrThrow<T>(name, input));
    }
    values = inputs_parsed;
  };
}

/**
 * Add an option to the CLI that accepts an enum value.
 * @param app The CLI app to add the option to
 * @param name The name of the option
 * @param value The value to store the option in
 * @param description The description of the option
 */
template <class T>
CLI::Option* AddOptionalEnumOption(AppPtr app, const std::string& name, T& value, const std::string& description) {
  std::string desc = fmt::format("value in {}", FormatEnumNames<T>("{", "}"));
  return app->add_option_function(name, EnumOptionFunction(name, value), description)
      ->run_callback_for_default()
      ->transform([](const std::string& value) -> std::string { return value; }, desc);
}

/**
 * Add an option to the CLI that accepts an enum value.
 * @param app The CLI app to add the option to
 * @param name The name of the option
 * @param value The value to store the option in
 * @param description The description of the option
 */
template <class T>
CLI::Option* AddOptionalEnumOption(AppPtr app,
                                   const std::string& name,
                                   std::vector<T>& value,
                                   const std::string& description) {
  std::string desc = fmt::format("value in {}", FormatEnumNames<T>("{", "}"));
  return app->add_option_function(name, EnumOptionFunction(name, value), description)
      ->delimiter(',')
      ->run_callback_for_default()
      ->transform([](const std::string& value) -> std::string { return value; }, desc);
}

/**
 * Add an option to the CLI that accepts an enum value.
 * @param app The CLI app to add the option to
 * @param name The name of the option
 * @param value The value to store the option in
 * @param description The description of the option
 * @param default_value The default value of the option
 */
template <class T>
CLI::Option* AddEnumOption(
    AppPtr app, const std::string& name, T& value, const std::string& description, const T& default_value) {
  return AddOptionalEnumOption(app, name, value, description)->default_val(enum_util::FormatEnumName(default_value));
}

/**
 * Add an option to the CLI that accepts an enum value.
 * @param app The CLI app to add the option to
 * @param name The name of the option
 * @param value The value to store the option in
 * @param description The description of the option
 * @param default_value The default value of the option
 */
template <class T>
CLI::Option* AddEnumOption(AppPtr app,
                           const std::string& name,
                           std::vector<T>& value,
                           const std::string& description,
                           const std::vector<T>& default_value) {
  return AddOptionalEnumOption(app, name, value, description)->default_val(FormatEnumName(default_value));
}

}  // namespace xoos::cli
