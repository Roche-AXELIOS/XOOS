#include "xoos/util/sequence-functions.h"

#include <ranges>

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

}  // namespace xoos::sequence
