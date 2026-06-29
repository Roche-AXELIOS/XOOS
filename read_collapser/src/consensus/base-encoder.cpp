#include "consensus/base-encoder.h"

namespace xoos::read_collapser {

BaseIndex ToIndex(const char base) {
  using enum BaseIndex;
  switch (base) {
    case kBaseA:
      return kIndexA;
    case kBaseC:
      return kIndexC;
    case kBaseG:
      return kIndexG;
    case kBaseT:
      return kIndexT;
    case kBaseGap:
    default:
      // Return kGap for invalid base
      return kIndexGap;
  }
}

size_t ToSizeTIndex(const char base) {
  return static_cast<size_t>(ToIndex(base));
}

char ToBase(const BaseIndex index) {
  using enum BaseIndex;
  switch (index) {
    case kIndexA:
      return kBaseA;
    case kIndexC:
      return kBaseC;
    case kIndexG:
      return kBaseG;
    case kIndexT:
      return kBaseT;
    case kIndexGap:
    default:
      // Return kBaseGap for invalid index
      return kBaseGap;
  }
}

}  // namespace xoos::read_collapser
