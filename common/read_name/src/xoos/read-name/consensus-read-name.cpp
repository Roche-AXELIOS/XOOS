#include "xoos/read-name/consensus-read-name.h"

#include <charconv>
#include <cstddef>
#include <exception>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/util/string-functions.h>

namespace xoos::read_name {

namespace {
// Index of the average-family-size field within the hyphen-delimited fields.
constexpr size_t kFamilySizeFieldIndex = 5;
}  // namespace

bool IsConsensusReadName(const std::string_view& qname) {
  const auto parts = string::SplitStringView(qname, '-');
  if (parts.size() != kConsensusReadNameFieldCount) {
    return false;
  }
  // Cluster id must be `{super_region_id}:{count}` with both sides numeric.
  const std::string_view cluster_id = parts.front();
  const size_t colon_pos = cluster_id.find(':');
  if (colon_pos == std::string_view::npos || colon_pos == 0 || colon_pos + 1 >= cluster_id.size() ||
      cluster_id.find(':', colon_pos + 1) != std::string_view::npos) {
    return false;
  }
  if (!string::OnlyDigits(cluster_id.substr(0, colon_pos)) || !string::OnlyDigits(cluster_id.substr(colon_pos + 1))) {
    return false;
  }
  // The last five fields (all but the cluster id) must be numeric.
  for (size_t i = 1; i < parts.size(); ++i) {
    if (!string::OnlyDigits(parts[i])) {
      return false;
    }
  }
  return true;
}

ConsensusReadInfo ParseConsensusReadName(const std::string_view& qname) {
  if (!IsConsensusReadName(qname)) {
    throw error::Error("'{}' does not follow the consensus read name format", qname);
  }
  const auto parts = string::SplitStringView(qname, '-');
  const auto parse_u32 = [&qname](const std::string_view& field) {
    u32 value = 0;
    const auto result = std::from_chars(field.data(), field.data() + field.size(), value, 10);
    if (result.ec != std::errc() || result.ptr != field.data() + field.size()) {
      throw error::Error("'{}' does not contain a valid consensus read name field '{}'", qname, field);
    }
    return value;
  };
  return {
      .forward_partial = parse_u32(parts[1]),
      .reverse_partial = parse_u32(parts[2]),
      .forward_full = parse_u32(parts[3]),
      .reverse_full = parse_u32(parts[4]),
      .avg_family_size = parse_u32(parts[kFamilySizeFieldIndex]),
  };
}

std::optional<ConsensusReadInfo> ParseConsensusReadNameNoExcept(const std::string_view& qname) noexcept {
  try {
    return ParseConsensusReadName(qname);
  } catch (const std::exception&) {
    return std::nullopt;
  }
}

std::string FormatConsensusReadName(const u32 super_region_id, const u32 count, const ConsensusReadInfo& info) {
  return fmt::format("{}:{}-{}-{}-{}-{}-{}",
                     super_region_id,
                     count,
                     info.forward_partial,
                     info.reverse_partial,
                     info.forward_full,
                     info.reverse_full,
                     info.avg_family_size);
}

}  // namespace xoos::read_name
