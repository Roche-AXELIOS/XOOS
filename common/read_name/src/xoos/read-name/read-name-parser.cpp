#include "xoos/read-name/read-name-parser.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <exception>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>

#include <xoos/error/error.h>
#include <xoos/types/int.h>

#include "xoos/read-name/umi.h"

namespace xoos::read_name {

bool IsLegacyReadNameFormat(const std::string_view& qname) {
  const size_t last_pipe_pos = qname.rfind('|');
  if (last_pipe_pos == std::string_view::npos) {
    throw error::Error("{} does not contain '|'", qname);
  }
  if (last_pipe_pos == 0) {
    throw error::Error("{} does not contain valid read name format", qname);
  }

  // Check if there's a second pipe before the last one
  const size_t second_last_pipe_pos = qname.rfind('|', last_pipe_pos - 1);

  // If there's no second pipe, it could be the new format (uses ':' as delimiter)
  // If there are 2+ pipes, it's the legacy format (uses '|' as delimiter)
  return (second_last_pipe_pos != std::string_view::npos);
}

bool HasUmiFields(const std::string_view& qname) noexcept {
  const size_t last_pipe_pos = qname.rfind('|');
  if (last_pipe_pos == std::string_view::npos || last_pipe_pos == 0) {
    return false;
  }
  // Legacy format (2+ pipes) always carries UMI fields.
  const size_t second_last_pipe_pos = qname.rfind('|', last_pipe_pos - 1);
  if (second_last_pipe_pos != std::string_view::npos) {
    return true;
  }
  // Newer format: UMI fields are present when the section after the last '|' contains a ':'
  // (`<sid>:<umi5p>:<umi3p>`). A bare `<sid>` indicates no UMI fields.
  const std::string_view after_last_pipe = qname.substr(last_pipe_pos + 1);
  return after_last_pipe.find(':') != std::string_view::npos;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity) - sequential read-name format parsing
AdapterInfo ParseReadName(const std::string_view& qname, const bool parse_umi) {
  // Detect format by counting pipes
  // Legacy format: <read_info>|<sid>|<umi5p>|<umi3p>
  // Newer format: <read_info>:<bitfield>|<sid>:<umi5p>:<umi3p> (`:<umi5p>:<umi3p>` is optional if UMIs are not present)
  // Newer format no UMIs: <read_info>|<sid>
  const size_t last_pipe_pos = qname.rfind('|');
  if (last_pipe_pos == std::string_view::npos) {
    throw error::Error("{} does not contain '|'", qname);
  }
  const bool is_legacy_format = IsLegacyReadNameFormat(qname);
  const char umi_delimiter = is_legacy_format ? '|' : ':';
  Umi umi5p = kNoUmi;
  Umi umi3p = kNoUmi;
  if (parse_umi) {
    // Find the last occurrence of the delimiter within the UMI section
    const size_t last_umi_delim_pos = qname.rfind(umi_delimiter);

    // If the delimiter is not found, it means there are no UMIs in the input
    if (last_umi_delim_pos == std::string_view::npos || last_umi_delim_pos == 0) {
      throw error::Error("'{}' does not contain expected UMI delimiters '{}'", qname, umi_delimiter);
    }

    // Find the second-to-last occurrence of the delimiter within the UMI section
    const size_t second_last_umi_delim_pos = qname.rfind(umi_delimiter, last_umi_delim_pos - 1);
    if (second_last_umi_delim_pos == std::string_view::npos) {
      throw error::Error("'{}' does not contain two UMI delimiters '{}'", qname, umi_delimiter);
    }
    umi5p =
        std::string(qname.substr(second_last_umi_delim_pos + 1, last_umi_delim_pos - second_last_umi_delim_pos - 1));
    umi3p = std::string(qname.substr(last_umi_delim_pos + 1));
    // Treat empty string UMIs as no UMIs.
    umi5p = umi5p->empty() ? kNoUmi : umi5p;
    umi3p = umi3p->empty() ? kNoUmi : umi3p;
  }

  // Extract SID and UMI presence info from the read name based on the adapter flag bits if it's the newer format
  bool has_5p_sid = true;
  bool has_3p_sid = true;
  bool has_5p_umi = false;
  bool has_3p_umi = false;
  // Only parse the adapter flag if the format is the new format.
  if (!is_legacy_format) {
    // In the new format, the adapter flag is encoded in the field immediately preceding the last pipe (`|`)
    // and is separated from the read info by a colon (`:`). The adapter flag is a bitfield where:
    // - bit 0 (1): presence of 3' UMI
    // - bit 1 (2): presence of 5' UMI
    // - bit 2 (4): presence of 3' SID
    // - bit 3 (8): presence of 5' SID
    const size_t adapter_flag_delim_pos = qname.rfind(':', last_pipe_pos - 1);
    if (adapter_flag_delim_pos == std::string_view::npos) {
      throw error::Error("'{}' does not contain expected adapter flag delimiter ':'", qname);
    }
    const std::string_view adapter_flag_str =
        qname.substr(adapter_flag_delim_pos + 1, last_pipe_pos - adapter_flag_delim_pos - 1);
    if (std::ranges::find_if_not(adapter_flag_str, ::isdigit) != adapter_flag_str.end()) {
      throw error::Error("'{}' contains non-numeric characters in the adapter flag field", qname);
    }
    u8 parsed_adapter_flag_value = 0;
    const auto parsing_result = std::from_chars(
        adapter_flag_str.data(), adapter_flag_str.data() + adapter_flag_str.size(), parsed_adapter_flag_value, 10);
    if (parsing_result.ec != std::errc()) {
      throw error::Error("'{}' does not contain a valid adapter bit field", qname);
    }
    has_5p_sid = (parsed_adapter_flag_value & kAdapterFlagHas5pSid) != 0;
    has_3p_sid = (parsed_adapter_flag_value & kAdapterFlagHas3pSid) != 0;
    has_5p_umi = (parsed_adapter_flag_value & kAdapterFlagHas5pUmi) != 0;
    has_3p_umi = (parsed_adapter_flag_value & kAdapterFlagHas3pUmi) != 0;
  }
  if (!is_legacy_format && parse_umi && ((has_5p_umi != (umi5p != kNoUmi)) || (has_3p_umi != (umi3p != kNoUmi)))) {
    throw error::Error("'{}' has inconsistent adapter flag and UMI presence", qname);
  }
  return {
      .has_5p_sid = has_5p_sid,
      .has_3p_sid = has_3p_sid,
      .has_5p_umi = is_legacy_format ? umi5p != kNoUmi : has_5p_umi,
      .has_3p_umi = is_legacy_format ? umi3p != kNoUmi : has_3p_umi,
      .umi5p = umi5p != kNoUmi ? umi5p : std::nullopt,
      .umi3p = umi3p != kNoUmi ? umi3p : std::nullopt,
  };
}

AdapterInfo ParseReadNameNoExcept(const std::string_view& qname, const bool parse_umi) noexcept {
  try {
    return ParseReadName(qname, parse_umi);
  } catch (const std::exception&) {
    return {
        .has_5p_sid = true,
        .has_3p_sid = true,
        .has_5p_umi = false,
        .has_3p_umi = false,
        .umi5p = std::nullopt,
        .umi3p = std::nullopt,
    };
  }
}

}  // namespace xoos::read_name
