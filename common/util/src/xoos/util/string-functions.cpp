#include "xoos/util/string-functions.h"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <locale>
#include <string>
#include <string_view>
#include <vector>

#include <xoos/error/error.h>

namespace xoos::string {

std::vector<std::string> Split(const std::string& s, const std::string& delimiter) {
  size_t pos_start = 0;
  size_t pos_end = 0;
  const size_t delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;
  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr(pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res.push_back(token);
  }
  res.push_back(s.substr(pos_start));
  return res;
}

std::vector<std::string_view> SplitStringView(const std::string_view s, const char delimiter) {
  std::vector<std::string_view> res;
  size_t start = 0;
  while (true) {
    const size_t pos = s.find(delimiter, start);
    if (pos == std::string_view::npos) {
      res.emplace_back(s.substr(start));
      break;
    }
    res.emplace_back(s.substr(start, pos - start));
    start = pos + 1;
  }
  return res;
}

bool OnlyDigits(const std::string_view s) {
  return !s.empty() && std::ranges::all_of(s, [](const char c) { return std::isdigit(c, std::locale::classic()); });
}

std::string BoolToString(const bool value) {
  return value ? "true" : "false";
}

bool StringToBool(const std::string& value) {
  if (value == "1" || value == "TRUE" || value == "true") {
    return true;
  }
  if (value == "0" || value == "FALSE" || value == "false") {
    return false;
  }
  throw error::Error("Invalid boolean string: {}", value);
}

std::string Join(const std::vector<std::string>& input_vector, const std::string& delimiter) {
  std::string result;
  uint32_t counter = 0;
  for (const auto& element : input_vector) {
    if (counter == input_vector.size() - 1) {
      result += element;
    } else {
      result += element + delimiter;
    }
    ++counter;
  }
  return result;
}

}  // namespace xoos::string
