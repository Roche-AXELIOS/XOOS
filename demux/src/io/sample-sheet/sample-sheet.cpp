#include "sample-sheet.h"

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include <csv.hpp>

namespace xoos::demux {
std::unordered_map<std::string, std::string> Read(const fs::path& path) {
  constexpr auto kSampleNameColumnHeader = "sample_name";
  constexpr auto kSampleSidColumnHeader = "sample_sid";

  auto sid2sample_name = std::unordered_map<std::string, std::string>{};
  auto sample_name2sid = std::unordered_map<std::string, std::string>{};
  auto reader = csv::CSVReader(path.string());
  for (auto& row : reader) {
    auto sample_name = row[kSampleNameColumnHeader].get<std::string>();
    auto sample_sid = row[kSampleSidColumnHeader].get<std::string>();

    auto it = sid2sample_name.find(sample_sid);
    if (it != sid2sample_name.end()) {
      throw error::Error("Found duplicate SID in sample sheet for '{}' and '{}'", sample_name, it->second);
    }
    sid2sample_name[sample_sid] = sample_name;

    auto it_name = sample_name2sid.find(sample_name);
    if (it_name != sample_name2sid.end()) {
      throw error::Error("Found duplicate sample name in sample sheet for '{}' and '{}'", sample_sid, it_name->second);
    }
    sample_name2sid[sample_name] = sample_sid;
  }
  return sid2sample_name;
}
}  // namespace xoos::demux
