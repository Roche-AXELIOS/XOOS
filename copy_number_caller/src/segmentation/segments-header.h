#pragma once
#include <string>
#include <string_view>

namespace xoos::cnc {
// Custom metadata keys of the SEG header, in the XOOS canonical PascalCase spelling. SEG files written by releases
// predating that convention spell them in UPPER_SNAKE_CASE; ReadHeaderFromSegments accepts both.
inline constexpr std::string_view kSegmentsHeaderProgramName = "ProgramName";
inline constexpr std::string_view kSegmentsHeaderVersion = "Version";
inline constexpr std::string_view kSegmentsHeaderReferenceFile = "ReferenceFile";
inline constexpr std::string_view kSegmentsHeaderCommandLine = "CommandLine";

struct SegmentsHeader {
  std::string program_name;
  std::string version;
  std::string reference_file;
  std::string command_line;
};
}  // namespace xoos::cnc
