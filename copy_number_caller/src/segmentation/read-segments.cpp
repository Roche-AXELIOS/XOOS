#include "segmentation/read-segments.h"

#include <algorithm>
#include <locale>
#include <string>
#include <string_view>
#include <utility>

#include <fmt/core.h>

#include <csv.hpp>

#include <xoos/enum/enum-util.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/vec.h>
#include <xoos/util/string-functions.h>

#include "likelihood/likelihood.h"
#include "segmentation/genomic-segment-column-names.h"
#include "segmentation/segments-header.h"

namespace xoos::cnc {

/**
 * @brief Extract the value of a CLI option from a serialized command line string.
 *
 * This parses the already-rendered command line stored in ##RocheCommandLine metadata.
 * CLI11 cannot be used here because we are operating on a flat string, not a live CLI11 app.
 * The format is always "--option value" with space separation (produced by RenderCli).
 */
namespace {
std::string ExtractCliOption(const std::string& command_line, const std::string& option_name) {
  auto tokens = xoos::string::Split(command_line, " ");
  auto it = std::ranges::find(tokens, option_name);
  if (it == tokens.end()) {
    return "";
  }
  ++it;
  if (it == tokens.end()) {
    return "";
  }
  return *it;
}

/**
 * @brief Fold a SEG metadata key into a spelling-insensitive form for comparison.
 *
 * Uppercasing and dropping underscores collapses the canonical PascalCase key (ProgramName) and the
 * UPPER_SNAKE_CASE key written by releases predating the metadata convention (PROGRAM_NAME) onto the same
 * token, so both spellings resolve to the same header field.
 */
std::string NormalizeHeaderKey(const std::string_view key) {
  const std::locale& classic = std::locale::classic();
  std::string normalized;
  normalized.reserve(key.size());
  for (const char c : key) {
    if (c != '_') {
      normalized.push_back(std::toupper(c, classic));
    }
  }
  return normalized;
}

vec<std::string> GetColumnNames(SegmentType mode) {
  using enum SegmentType;
  switch (mode) {
    case kLogROnly:
      return kGenomicSegColsRequiredLogRSegments;
    case kBaf:
      return kGenomicSegColsRequiredBafSegments;
    case kGermlineLikelihood:
      return kGenomicSegColsRequiredGermline;
    case kSomaticWithBafLikelihood:
      return kGenomicSegColsRequiredSomaticAlleleSpecific;
    case kSomaticNoBafLikelihood:
      return kGenomicSegColsRequiredSomaticNoAlleleSpecific;
    case kSomaticTumorOnlyWithBafLikelihood:
      return kGenomicSegColsRequiredSomaticTumorOnlyAlleleSpecific;
    case kSeed:
      return kGenomicSegColsRequiredSeed;
    default:
      Logging::Error("invalid SegmentType {}", enum_util::FormatEnumName(mode));
      throw std::runtime_error("invalid SegmentType");
  }
}
}  // namespace

/** @brief
 *
 * we have to read a string here so that the parser can "magically" determine the format and know when to skip comments
 */
vec<GenomicSegment> ReadSegments(const fs::path& fname, SegmentType segment_type) {
  vec<GenomicSegment> ret;
  csv::CSVReader reader(fname.string());
  auto col_names = GetColumnNames(segment_type);
  s32 row_num = 0;
  for (const csv::CSVRow& row : reader) {
    // check if all the required columns are present in this row
    // csv-parser 5.x returns col names by const reference (2.3.0 returned by value); bind to avoid a copy.
    const auto& row_col_names = row.get_col_names();  // NOSONAR(M23_225) — binds a const reference; no allocation
    GenomicSegment seg;
    for (const std::string& col : col_names) {
      if (std::ranges::find(row_col_names, col) == row_col_names.end()) {
        Logging::Error("missing column {} in file {} at row number {}", col, fname.string(), row_num);
        throw std::runtime_error(fmt::format("missing column {}", col));
      }
      csv::CSVField val = row[col];
      seg.SetFieldOrThrowError(col, val.get<std::string>());
      if (segment_type == SegmentType::kGermlineLikelihood) {
        seg.purity = kGermlinePurity;
        seg.ploidy = kDiploidPloidy;
      }
    }
    if (seg.end <= seg.start) {
      Logging::Error("SEG file error: end position must be greater than start position in file {} at row number {}",
                     fname.string(),
                     row_num);
      throw std::runtime_error("invalid segment");
    }
    ret.emplace_back(seg);
    row_num += 1;
  }
  Logging::Debug("read {} rows", row_num);
  return ret;
}

SegmentsHeader ReadHeaderFromSegments(const fs::path& fname) {
  const std::string program_name_key = NormalizeHeaderKey(kSegmentsHeaderProgramName);
  const std::string version_key = NormalizeHeaderKey(kSegmentsHeaderVersion);
  const std::string reference_file_key = NormalizeHeaderKey(kSegmentsHeaderReferenceFile);
  const std::string command_line_key = NormalizeHeaderKey(kSegmentsHeaderCommandLine);

  SegmentsHeader ret;
  std::ifstream ifs(fname);
  std::string line;
  while (std::getline(ifs, line)) {
    if (!line.starts_with('#')) {
      break;
    }
    // New format: ##RocheCommandLine=<ID=...,Version="...",CommandLine="...">
    auto parsed = io::ParseCommandInfo(line);
    if (parsed.has_value()) {
      ret.program_name = parsed->name;
      ret.version = parsed->version;
      ret.command_line = parsed->command_line;
      if (ret.reference_file.empty()) {
        ret.reference_file = ExtractCliOption(parsed->command_line, "--reference");
      }
      continue;
    }
    // Custom metadata line: ##Key=Value. Releases predating the metadata convention wrote a single '#', an
    // UPPER_SNAKE_CASE key, and either '=' or ':' as the delimiter, so all of those spellings are accepted.
    const size_t key_begin = line.find_first_not_of('#');
    if (key_begin == std::string::npos) {
      continue;
    }
    const size_t delim = line.find_first_of("=:", key_begin);
    if (delim == std::string::npos) {
      continue;
    }
    const std::string key = NormalizeHeaderKey(std::string_view{line}.substr(key_begin, delim - key_begin));
    std::string value = line.substr(delim + 1);
    if (key == program_name_key) {
      ret.program_name = std::move(value);
    } else if (key == version_key) {
      ret.version = std::move(value);
    } else if (key == reference_file_key) {
      ret.reference_file = std::move(value);
    } else if (key == command_line_key) {
      ret.command_line = std::move(value);
    }
  }
  return ret;
}
}  // namespace xoos::cnc
