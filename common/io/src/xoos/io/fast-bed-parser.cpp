#include "xoos/io/fast-bed-parser.h"

#include <cctype>
#include <charconv>
#include <cstddef>
#include <locale>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <htslib/hts.h>
#include <htslib/kseq.h>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/types/int.h>

namespace xoos::io {

namespace {

constexpr size_t kMinimumBedColumns = 3;
constexpr size_t kMinimumBedColumnsWithData = 4;

bool IsBedWhitespace(const char character) {
  return std::isspace(character, std::locale::classic());
}

std::string_view TrimBedField(std::string_view field) {
  const size_t first = field.find_first_not_of(" \t\r\n\v\f");
  if (first == std::string_view::npos) {
    return {};
  }
  const size_t last = field.find_last_not_of(" \t\r\n\v\f");
  return field.substr(first, last - first + 1);
}

struct BedField {
  std::string_view value;
  size_t field_end;
};

std::optional<BedField> GetNextBedField(const std::string_view line, size_t field_start) {
  while (field_start < line.size() && IsBedWhitespace(line[field_start])) {
    ++field_start;
  }
  if (field_start == line.size()) {
    return std::nullopt;
  }
  size_t field_end = field_start;
  while (field_end < line.size() && !IsBedWhitespace(line[field_end])) {
    ++field_end;
  }
  return BedField{.value = line.substr(field_start, field_end - field_start), .field_end = field_end};
}

bool IsBedComment(const std::string_view field) {
  return field.starts_with('#');
}

bool IsBedDirective(const std::string_view field) {
  return field == "browser" || field == "track";
}

bool ParseBedCoordinate(const std::string_view text, s64& coordinate) {
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  // std::from_chars uses raw pointers to the character array of the string
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), coordinate);
  return error == std::errc{} && end == text.data() + text.size();
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

struct BedCoordinates {
  s64 start;
  s64 end;
};

std::optional<BedCoordinates> ParseBedCoordinates(const BedField& start_field,
                                                  const BedField& end_field,
                                                  const std::string_view chromosome) {
  s64 start = 0;
  s64 end = 0;
  if (!ParseBedCoordinate(start_field.value, start) || !ParseBedCoordinate(end_field.value, end)) {
    if (IsBedDirective(chromosome)) {
      return std::nullopt;
    }
    throw error::Error("BED file must have numeric start and end columns");
  }
  if (start < 0 || end < 0) {
    if (IsBedDirective(chromosome)) {
      return std::nullopt;
    }
    throw error::Error("BED file must have non-negative start and end columns");
  }
  if (start > end) {
    if (IsBedDirective(chromosome)) {
      return std::nullopt;
    }
    throw error::Error("BED file must have start <= end");
  }
  return BedCoordinates{.start = start, .end = end};
}

std::vector<std::string> ParseTabDelimitedBedExtraColumns(const std::string_view line,
                                                          const size_t first_column_start) {
  std::vector<std::string> data;
  size_t field_start = first_column_start;
  while (true) {
    const size_t field_end = line.find('\t', field_start);
    const size_t field_length = field_end == std::string_view::npos ? std::string_view::npos : field_end - field_start;
    data.emplace_back(TrimBedField(line.substr(field_start, field_length)));
    if (field_end == std::string_view::npos) {
      return data;
    }
    field_start = field_end + 1;
  }
}

std::vector<std::string> ParseWhitespaceDelimitedBedExtraColumns(const std::string_view line, size_t field_start) {
  std::vector<std::string> data;
  while (const std::optional<BedField> field = GetNextBedField(line, field_start)) {
    data.emplace_back(field->value);
    field_start = field->field_end;
  }
  return data;
}

std::vector<std::string> ParseBedExtraColumns(const std::string_view line, const size_t third_field_end) {
  if (third_field_end < line.size() && line[third_field_end] == '\t') {
    return ParseTabDelimitedBedExtraColumns(line, third_field_end + 1);
  }
  return ParseWhitespaceDelimitedBedExtraColumns(line, third_field_end);
}

[[noreturn]] void ThrowInsufficientColumns(const bool expect_extra_columns) {
  if (expect_extra_columns) {
    throw error::Error(
        "BED file must have at least {} columns in order to additionally parse data from the fourth column onward",
        kMinimumBedColumnsWithData);
  }
  throw error::Error("BED file must have at least {} columns", kMinimumBedColumns);
}

}  // namespace

std::optional<BedRegion> ParseBedRowFast(const std::string_view line, const bool expect_extra_columns) {
  const std::optional<BedField> chromosome_field = GetNextBedField(line, 0);
  if (!chromosome_field.has_value()) {
    ThrowInsufficientColumns(expect_extra_columns);
  }
  const std::string_view chromosome = chromosome_field->value;
  if (IsBedComment(chromosome)) {
    return std::nullopt;
  }
  const std::optional<BedField> start_field = GetNextBedField(line, chromosome_field->field_end);
  const std::optional<BedField> end_field =
      start_field.has_value() ? GetNextBedField(line, start_field->field_end) : std::nullopt;
  if (!start_field.has_value() || !end_field.has_value()) {
    if (IsBedDirective(chromosome)) {
      return std::nullopt;
    }
    ThrowInsufficientColumns(expect_extra_columns);
  }

  const std::optional<BedCoordinates> coordinates = ParseBedCoordinates(*start_field, *end_field, chromosome);
  if (!coordinates.has_value()) {
    return std::nullopt;
  }
  const std::optional<BedField> first_extra_field = GetNextBedField(line, end_field->field_end);
  const bool has_tab_delimited_extra_column = end_field->field_end < line.size() && line[end_field->field_end] == '\t';
  if (expect_extra_columns && !first_extra_field.has_value() && !has_tab_delimited_extra_column) {
    ThrowInsufficientColumns(true);
  }

  std::vector<std::string> data;
  if (first_extra_field.has_value() || has_tab_delimited_extra_column) {
    data = ParseBedExtraColumns(line, end_field->field_end);
  }
  return BedRegion{std::string(chromosome), coordinates->start, coordinates->end, std::move(data)};
}

}  // namespace xoos::io
