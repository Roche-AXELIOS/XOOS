#include "xoos/util/sequence-functions.h"

#include <algorithm>
#include <ranges>
#include <string>
#include <string_view>

#include <xoos/types/int.h>

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
      str.crbegin(), str.crend(), reverse_complement.begin(), [](const u8 c) { return kBaseToComplement.at(c); });
  return reverse_complement;
}

bool IsACGT(const char c) {
  return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

bool IsAmbiguousBase(const char c) {  // NOSONAR cpp:S1541 - switch statement used for efficient character matching
  switch (c) {
    case 'R':
    case 'r':
    case 'Y':
    case 'y':
    case 'S':
    case 's':
    case 'W':
    case 'w':
    case 'K':
    case 'k':
    case 'M':
    case 'm':
    case 'B':
    case 'b':
    case 'D':
    case 'd':
    case 'H':
    case 'h':
    case 'V':
    case 'v':
      return true;
    default:
      return false;
  }
}

bool AmbiguousBaseMatchesQuery(  // NOSONAR cpp:S1541 - switch statement used for efficient character matching
    const char ref_upper,
    const char query_upper) {
  switch (ref_upper) {
    case 'R':
      return query_upper == 'A' || query_upper == 'G';
    case 'Y':
      return query_upper == 'C' || query_upper == 'T';
    case 'S':
      return query_upper == 'G' || query_upper == 'C';
    case 'W':
      return query_upper == 'A' || query_upper == 'T';
    case 'K':
      return query_upper == 'G' || query_upper == 'T';
    case 'M':
      return query_upper == 'A' || query_upper == 'C';
    case 'B':
      return query_upper == 'C' || query_upper == 'G' || query_upper == 'T';
    case 'D':
      return query_upper == 'A' || query_upper == 'G' || query_upper == 'T';
    case 'H':
      return query_upper == 'A' || query_upper == 'C' || query_upper == 'T';
    case 'V':
      return query_upper == 'A' || query_upper == 'C' || query_upper == 'G';
    default:
      return false;
  }
}

}  // namespace xoos::sequence
