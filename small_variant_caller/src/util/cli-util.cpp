#include "cli-util.h"

#include <functional>

#include <fmt/format.h>

#include <xoos/enum/enum-util.h>

#include "core/cli-option-names.h"
#include "log-util.h"

namespace xoos::svc {

CLI::Option* AddWarnAsErrorOption(CLI::App* const app) {
  return app->add_flag(cli_opt_name::kWarnAsError, warn_as_error, "Treat warn messages as errors");
}

CLI::Option* AddExplicitLibraryTypeOption(CLI::App* const app,
                                          LibraryType& library_type,
                                          const std::string& description) {
  constexpr std::string_view kExplicitLibraryTypeValues = "{sbx-s, sbx-d, sbx-fast}";
  const auto validation_message = fmt::format("value in {}", kExplicitLibraryTypeValues);
  const std::function<void(const std::string&)> parse_library_type =
      [&library_type, kExplicitLibraryTypeValues](const std::string& input) {
        const auto parsed = enum_util::ParseEnumName<LibraryType>(input);
        if (!parsed.has_value() || *parsed == LibraryType::kAutoDetect) {
          throw CLI::ValidationError(fmt::format(
              "{}: Check {} value in {} FAILED", cli_opt_name::kLibraryType, input, kExplicitLibraryTypeValues));
        }
        library_type = *parsed;
      };
  return app->add_option_function(cli_opt_name::kLibraryType, parse_library_type, description)
      ->run_callback_for_default()
      ->transform([](const std::string& value) -> std::string { return value; }, validation_message);
}

void SetRenderedDefault(CLI::App* const sub, const std::string_view option_name, const std::string& value) {
  auto* const opt = sub->get_option_no_throw(std::string{option_name});
  if (opt != nullptr && opt->empty()) {
    opt->default_str(value);
  }
}

void ReplaceOrAppendCommandLineOption(std::string& command_line,
                                      const std::string_view option_name,
                                      const std::string& value) {
  const std::string option{option_name};

  // Match the flag as a whole token so a shorter flag does not match the prefix of a longer one
  // (e.g. "--sample-type" must not match inside "--sample-type-max-af"). A token boundary requires
  // the flag to be preceded by the start of the string or whitespace and followed by whitespace.
  auto option_pos = std::string::npos;
  for (auto pos = command_line.find(option); pos != std::string::npos; pos = command_line.find(option, pos + 1)) {
    const bool starts_token = pos == 0 || command_line[pos - 1] == ' ';
    const auto after = pos + option.size();
    const bool ends_token = after >= command_line.size() || command_line[after] == ' ';
    if (starts_token && ends_token) {
      option_pos = pos;
      break;
    }
  }

  if (option_pos == std::string::npos) {
    command_line += fmt::format(" {} {}", option, value);
    return;
  }

  const auto value_start = command_line.find_first_not_of(' ', option_pos + option.size());
  if (value_start == std::string::npos) {
    command_line += fmt::format(" {}", value);
    return;
  }
  const auto value_end = command_line.find(' ', value_start);
  command_line.replace(
      value_start, value_end == std::string::npos ? std::string::npos : value_end - value_start, value);
}

}  // namespace xoos::svc
