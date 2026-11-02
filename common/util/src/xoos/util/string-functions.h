#pragma once

#include <algorithm>
#include <locale>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

namespace xoos::string {

// string splitting function missing in standard library
std::vector<std::string> Split(const std::string& s, const std::string& delimiter);

// Split a string_view on a single-character delimiter, returning views into the original string
// without allocating a copy of each field. The input must outlive the returned views. An empty input
// yields a single empty field, mirroring `Split`.
std::vector<std::string_view> SplitStringView(std::string_view s, char delimiter);

// Return true if the string is non-empty and contains only decimal digit characters.
bool OnlyDigits(std::string_view s);

// Convert a boolean to the standard lowercase value-representation spelling.
std::string BoolToString(bool value);

// Parse a tolerant historical boolean spelling used in existing XOOS text outputs.
bool StringToBool(const std::string& value);

// trim from start (in place)
static inline void LeftTrim(std::string& s) {
  s.erase(s.begin(), std::ranges::find_if(s, [](const char ch) { return !std::isspace(ch, std::locale::classic()); }));
}

// trim from end (in place)
static inline void RightTrim(std::string& s) {
  s.erase(std::ranges::find_if(std::ranges::reverse_view(s),
                               [](const char ch) { return !std::isspace(ch, std::locale::classic()); })
              .base(),
          s.end());
}

// trim from both ends (in place)
static inline void Trim(std::string& s) {
  RightTrim(s);
  LeftTrim(s);
}

// trim from start (copying)
static inline std::string LeftTrimCopy(std::string s) {
  LeftTrim(s);
  return s;
}

// trim from end (copying)
static inline std::string RightTrimCopy(std::string s) {
  RightTrim(s);
  return s;
}

// trim from both ends (copying)
static inline std::string TrimCopy(std::string s) {
  Trim(s);
  return s;
}

std::string Join(const std::vector<std::string>& input_vector, const std::string& delimiter);

static inline bool EndsWithOneOf(const std::string& name, const std::vector<std::string>& suffixes) {
  std::string ps{name};
  std::ranges::transform(ps, std::begin(ps), [](const char c) { return std::tolower(c, std::locale::classic()); });
  return std::ranges::any_of(suffixes, [&ps](const auto& suffix) { return ps.ends_with(suffix); });
}

}  // namespace xoos::string
