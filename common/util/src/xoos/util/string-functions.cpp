#include "xoos/util/string-functions.h"

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace xoos::string {

std::vector<std::string> Split(const std::string& s, const std::string& delimiter) {
  size_t pos_start = 0;
  size_t pos_end;
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

UmiType ParseUmis(const std::string& read_name, std::string& umi1, std::string& umi2) {
  auto fields = Split(read_name, "|");
  if (fields.size() < 3) {
    return kNone;
  }
  // GPU version
  // bunch:of:colon:separated:fields:SID|UMI1|UMI2
  // where UMI1 and UMI2 are integers >= 0
  // Nothing in the code actually interprets th UMIs, so we can just leave them as-is
  umi1 = *(fields.end() - 2);
  umi2 = *(fields.end() - 1);
  if (umi1.empty() || umi2.empty()) {
    throw std::runtime_error("Empty UMI detected");
  }
  if ((umi1.starts_with('-') || (isdigit(umi1[0]) != 0)) || (umi2.starts_with('-') || (isdigit(umi2[0]) != 0))) {
    // GPU version
    // parse them to ensure they are numbers
    if (umi1 != "*") {
      std::stoi(umi1);
    }
    if (umi2 != "*") {
      std::stoi(umi2);
    }
    return string::kGpu;
  }
  // demux_lut version
  // bunch:of:colon:separated:fields|SID1|SID2|UMI1|UMI2
  // where UMI1 and UMI2 are {ACGT}+|\*
  if (umi1.find_first_not_of("ACGT") != std::string::npos && umi1 != "*") {
    throw std::runtime_error("Invalid UMI detected");
  }
  if (umi2.find_first_not_of("ACGT") != std::string::npos && umi2 != "*") {
    throw std::runtime_error("Invalid UMI detected");
  }
  return kDemuxLut;
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
