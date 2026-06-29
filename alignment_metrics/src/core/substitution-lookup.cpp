#include "core/substitution-lookup.h"

#include <algorithm>
#include <array>
#include <tuple>

#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

/**
 * Returns the index for the substitution of base for ref.
 * A->C = 0
 * A->G = 1
 * A->T = 2
 * A->N = 3
 * C->A = 4
 * C->G = 5
 * C->T = 6
 * C->N = 7
 * G->A = 8
 * G->C = 9
 * G->T = 10
 * G->N = 11
 * T->A = 12
 * T->C = 13
 * T->G = 14
 * T->N = 15
 * N->A = 16
 * N->C = 17
 * N->G = 18
 * N->T = 19
 * unknown = 20
 */
u32 SubstitutionLookup::GetIndexForSubstitution(const char ref, const char base) {
  const auto base_index_iter = std::ranges::find(kBases, base);
  if (base_index_iter == kBases.end()) {
    return kUnknownSubstitutionIndex;
  }
  const auto ref_index_iter = std::ranges::find(kBases, ref);
  if (ref_index_iter == kBases.end()) {
    return kUnknownSubstitutionIndex;
  }
  const auto ref_index = ref_index_iter - kBases.begin();
  const auto base_index = base_index_iter - kBases.begin();
  return kSubstitutionLookupTable[ref_index][base_index];
}

// returns the reverse of GetIndexForSubstitution
// only used for the initial constructor
std::tuple<char, char> SubstitutionLookup::GetSubstitutionForIndex(const u32 index) {
  switch (index) {
    case 0:
      return {'A', 'C'};
    case 1:
      return {'A', 'G'};
    case 2:
      return {'A', 'T'};
    case 3:
      return {'A', 'N'};
    case 4:
      return {'C', 'A'};
    case 5:
      return {'C', 'G'};
    case 6:
      return {'C', 'T'};
    case 7:
      return {'C', 'N'};
    case 8:
      return {'G', 'A'};
    case 9:
      return {'G', 'C'};
    case 10:
      return {'G', 'T'};
    case 11:
      return {'G', 'N'};
    case 12:
      return {'T', 'A'};
    case 13:
      return {'T', 'C'};
    case 14:
      return {'T', 'G'};
    case 15:
      return {'T', 'N'};
    case 16:
      return {'N', 'A'};
    case 17:
      return {'N', 'C'};
    case 18:
      return {'N', 'G'};
    case 19:
      return {'N', 'T'};
    default:
      return {'N', 'N'};
  }
}

const std::array<char, 5> SubstitutionLookup::kBases{'A', 'C', 'G', 'T', 'N'};
const std::array<std::array<u32, 5>, 5> SubstitutionLookup::kSubstitutionLookupTable{{
    {3, 0, 1, 2, 3},                             // A -> {N, C, G, T, N}
    {4, 7, 5, 6, 7},                             // C -> {A, N, G, T, N}
    {8, 9, 11, 10, 11},                          // G -> {A, C, N, T, N}
    {12, 13, 14, 15, 15},                        // T -> {A, C, G, N, N}
    {16, 17, 18, 19, kUnknownSubstitutionIndex}  // N -> {A, C, G, T, N}
}};

}  // namespace xoos::alignment_metrics
