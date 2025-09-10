#pragma once

#include <algorithm>
#include <string>
#include <vector>

namespace xoos::string {

// string splitting function missing in standard library
std::vector<std::string> Split(const std::string& s, const std::string& delimiter);

// trim from start (in place)
static inline void LeftTrim(std::string& s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) { return std::isspace(ch) == 0; }));
}

// trim from end (in place)
static inline void RightTrim(std::string& s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return std::isspace(ch) == 0; }).base(), s.end());
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

enum UmiType {
  kGpu,
  kDemuxLut,
  kNone
};

UmiType ParseUmis(const std::string& read_name, std::string& umi1, std::string& umi2);

std::string Join(const std::vector<std::string>& input_vector, const std::string& delimiter);

static inline void FastUppercase(std::string& result) {
  // use a fast inlined version of toupper
  // this is not locale-aware, so it will only work for ASCII
  for (char& c : result) {
    c &= ~0x20;
  }
}

static inline bool EndsWithOneOf(const std::string& name, const std::vector<std::string>& suffixes) {
  std::string ps{name};
  std::transform(std::cbegin(ps), std::cend(ps), std::begin(ps), [](auto c) { return std::tolower(c); });
  return std::any_of(suffixes.cbegin(), suffixes.cend(), [&ps](const auto& suffix) { return ps.ends_with(suffix); });
}

}  // namespace xoos::string
