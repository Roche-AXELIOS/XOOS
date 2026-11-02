#include "xoos/io/metadata-util.h"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <xoos/log/logging.h>

namespace xoos::io {

namespace {

/// Minimum length of a quoted string: two surrounding double-quote characters.
constexpr size_t kMinQuotedLength = 2;
/// Advance by 2 when consuming a backslash-escaped character (the backslash + the escaped char).
constexpr size_t kEscapeSequenceLength = 2;

std::string EscapeMetadataValue(const std::string& value) {
  std::string escaped;
  escaped.reserve(value.size());
  for (const char c : value) {
    if (c == '\\' || c == '"') {
      escaped += '\\';
    }
    escaped += c;
  }
  return escaped;
}

/**
 * @brief Unescape a value that was escaped by EscapeMetadataValue.
 *
 * Reverses backslash-escaping of '"' and '\' characters.
 */
std::string UnescapeMetadataValue(const std::string& value) {
  std::string result;
  result.reserve(value.size());
  size_t i = 0;
  while (i < value.size()) {
    if (value[i] == '\\' && i + 1 < value.size()) {
      result += value[i + 1];
      i += kEscapeSequenceLength;
    } else {
      result += value[i];
      ++i;
    }
  }
  return result;
}

/**
 * @brief Strip surrounding double quotes from a value and unescape the interior.
 *
 * Returns the unescaped content between the quotes, or the original value if it is not quoted.
 */
std::string StripQuotes(const std::string& value) {
  if (value.size() >= kMinQuotedLength && value.front() == '"' && value.back() == '"') {
    return UnescapeMetadataValue(value.substr(1, value.size() - kMinQuotedLength));
  }
  return value;
}

/**
 * @brief Split the body of a RocheCommandLine tag by commas, respecting quoted values.
 *
 * Commas inside double-quoted strings (including escaped quotes) are not treated as delimiters.
 */
std::vector<std::string> SplitKeyValuePairs(const std::string& body) {
  std::vector<std::string> tokens;
  std::string current;
  bool in_quotes = false;
  size_t i = 0;
  while (i < body.size()) {
    if (body[i] == '\\' && in_quotes && i + 1 < body.size()) {
      current += body[i];
      current += body[i + 1];
      i += kEscapeSequenceLength;
    } else if (body[i] == '"') {
      in_quotes = !in_quotes;
      current += body[i];
      ++i;
    } else if (body[i] == ',' && !in_quotes) {
      tokens.push_back(current);
      current.clear();
      ++i;
    } else {
      current += body[i];
      ++i;
    }
  }
  if (!current.empty()) {
    tokens.push_back(current);
  }
  return tokens;
}

}  // namespace

std::string GetCommandInfo(const CommandLineInfo& info) {
  return fmt::format(R"(##RocheCommandLine=<ID={},Version="{}",CommandLine="{}">)",
                     info.name,
                     EscapeMetadataValue(info.version),
                     EscapeMetadataValue(info.command_line));
}

std::string GetSampleInfo(const ReadGroup& read_group) {
  std::vector<std::string> fields;
  const auto append_field = [&fields](const std::string& key, const std::string& value) {
    if (value.empty()) {
      return;
    }
    // Values are unquoted by default. Wrap in double quotes when the value contains a comma
    // (the field delimiter) or a double quote, so the value is unambiguous; embedded quotes
    // and backslashes are always escaped.
    const auto escaped = EscapeMetadataValue(value);
    if (value.contains(',') || value.contains('"')) {
      fields.push_back(fmt::format(R"({}="{}")", key, escaped));
    } else {
      fields.push_back(fmt::format("{}={}", key, escaped));
    }
  };
  append_field(kKeyName, read_group.sample_name);
  append_field(kKeyLibraryId, read_group.library);
  append_field(kKeyPlatform, read_group.platform);
  append_field(kKeyPlatformModel, read_group.platform_model);
  return fmt::format("{}{}>", kSampleLinePrefix, fmt::join(fields, ","));
}

std::optional<CommandLineInfo> ParseCommandInfo(const std::string_view line) {
  if (!line.starts_with(kRocheCommandLinePrefix)) {
    return std::nullopt;
  }
  if (line.empty() || line.back() != '>') {
    Logging::Warn("Malformed RocheCommandLine line: missing closing '>'");
    return std::nullopt;
  }

  const std::string body{line.substr(kRocheCommandLinePrefix.size(), line.size() - kRocheCommandLinePrefix.size() - 1)};
  const auto tokens = SplitKeyValuePairs(body);

  CommandLineInfo info;
  bool has_id = false;
  bool has_version = false;
  bool has_command_line = false;

  for (const auto& token : tokens) {
    const auto eq_pos = token.find('=');
    if (eq_pos == std::string::npos) {
      Logging::Warn("Malformed RocheCommandLine token (no '='): '{}'", token);
      return std::nullopt;
    }
    const auto key = token.substr(0, eq_pos);
    const auto value = token.substr(eq_pos + 1);

    if (key == kKeyId) {
      info.name = value;
      has_id = true;
    } else if (key == kKeyVersion) {
      info.version = StripQuotes(value);
      has_version = true;
    } else if (key == kKeyCommandLine) {
      info.command_line = StripQuotes(value);
      has_command_line = true;
    }
  }

  if (!has_id || !has_version || !has_command_line) {
    Logging::Warn("Malformed RocheCommandLine line: missing required field(s) (ID, Version, CommandLine)");
    return std::nullopt;
  }

  return info;
}

}  // namespace xoos::io
