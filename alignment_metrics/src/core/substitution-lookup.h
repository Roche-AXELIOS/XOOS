#pragma once

#include <array>
#include <tuple>

#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

// This is a lookup table for substitutions.
// This provides a more compact and efficient way to store substitutions.
// At the time of writing, this was more performant than using a map.
class SubstitutionLookup {
 public:
  static constexpr u32 kUnknownSubstitutionIndex = 20;

  // if the base is not A, C, G, T, N, we return this index
  static constexpr u32 kUnknownBaseIndex = 4;

  static constexpr u32 kMaxSubstitutionIndex = kUnknownSubstitutionIndex + 1;

  // Gets the internal index (encoding) for the given substitution.
  static u32 GetIndexForSubstitution(char ref, char base);

  // Gets the substitution (decodes) for the given index.
  static std::tuple<char, char> GetSubstitutionForIndex(u32 index);

 private:
  // internally defined constant maps A, C, G, T, N
  static const std::array<char, 5> kBases;
  // internally defined array of all possible substitutions for the given bases
  static const std::array<std::array<u32, 5>, 5> kSubstitutionLookupTable;
};

}  // namespace xoos::alignment_metrics
