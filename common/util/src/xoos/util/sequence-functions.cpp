#include "xoos/util/sequence-functions.h"

#include <algorithm>
#include <ranges>
#include <string>

namespace xoos::sequence {

// simple reverse-complement used for processing UMIs, which may be
// set to just "*" for partial reads
std::string ReverseComplement(const std::string& str) {
  std::string result;
  result.reserve(str.length());
  for (auto it : std::ranges::reverse_view(str)) {
    switch (it) {
      case 'A':
      case 'a':
        result.push_back('T');
        break;
      case 'C':
      case 'c':
        result.push_back('G');
        break;
      case 'G':
      case 'g':
        result.push_back('C');
        break;
      case 'T':
      case 't':
        result.push_back('A');
        break;
      case '*':
        result.push_back('*');
        break;
      default:
        result.push_back('N');
        break;
    }
  }
  return result;
}

std::string ReverseComplementCaseSensitive(const std::string_view& str) {
  auto reverse_complement = std::string(str.size(), '\0');
  std::ranges::transform(
      str.crbegin(), str.crend(), reverse_complement.begin(), [](const u8 c) { return kBaseToComplement[c]; });
  return reverse_complement;
}

bool IsACGT(const char c) {
  return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

}  // namespace xoos::sequence
