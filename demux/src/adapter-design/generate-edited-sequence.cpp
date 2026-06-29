#include "adapter-design/generate-edited-sequence.h"

#include <fmt/format.h>
#include <xoos/types/vec.h>

#include <unordered_set>

namespace xoos::demux {
void GenerateEditedSequence1(const std::string_view& alphabet, const std::string_view& sequence,
                             vec<std::string>& edited_sequences);

/**
 * Generate all possible edited sequences in 2 phases:
 *  1. For each edit distance, derive all edited sequences of edit distance 1 from the sequences of edit distance - 1
 *  2. Deduplicate any collisions in sequences, prioritizing least edit distance
 */
void GenerateEditedSequence(const std::string_view& alphabet, const std::string_view& sequence, u32 max_edit_distance,
                            const EditedSequenceCb& cb) {
  // edited_sequences contains a map from edit distance to all sequences at that edit distance,
  // this is required to ensure that if we have collisions between edited sequences we prioritize
  // the least edit distance
  vec<vec<std::string>> edited_sequences{max_edit_distance + 1};

  // Initialize for edit distance = 0, which is simply the initial sequence, this is the base case
  edited_sequences.at(0).emplace_back(sequence);

  // Start at edit distance = 1, and generate all edit distance 1 sequences from the
  // edit distance - 1 sequences
  for (u32 edit_distance = 1; edit_distance <= max_edit_distance; edit_distance++) {
    // For each sequence of edit distance - 1, we generate the sequences of
    // we generate sequences of edit distance 1 from those sequences,
    // this will give us sequences of edit distance
    for (const auto& edited_sequence : edited_sequences.at(edit_distance - 1)) {
      GenerateEditedSequence1(alphabet, edited_sequence, edited_sequences.at(edit_distance));
    }
  }

  // Track previously seen sequences to:
  //  1. Not produce the same edited sequence twice
  //  2. Ensure for duplicate edited sequences, the one with the least edit distance is preferred
  std::unordered_set<std::string> visited;

  // Start at edit distance 0 to prefer sequences with the least edit distance
  for (u32 edit_distance = 0; edit_distance < edited_sequences.size(); ++edit_distance) {
    for (const auto& edited_sequence : edited_sequences.at(edit_distance)) {
      if (visited.contains(edited_sequence)) {
        continue;
      }
      visited.insert(edited_sequence);
      cb(edit_distance, edited_sequence);
    }
  }
}

/*
 * Provided a sequence and an alphabet, this GenerateEditedSequence1 will
 * generate all sequences of edit distance 1. It does this by generating all possible
 * substitutions, deletions, and insertions at each location in the original sequence.
 */
void GenerateEditedSequence1(const std::string_view& alphabet, const std::string_view& sequence,
                             vec<std::string>& edited_sequences) {
  for (std::string::size_type i = 0; i < sequence.length() + 1; ++i) {
    std::string_view a = sequence.substr(0, i);
    std::string_view b = sequence.substr(i);

    // deletion
    edited_sequences.emplace_back(fmt::format("{}{}", a, b.empty() ? "" : b.substr(1)));
    for (auto c : alphabet) {
      // substitution
      edited_sequences.emplace_back(fmt::format("{}{}{}", a, c, b.empty() ? "" : b.substr(1)));
      // insert
      edited_sequences.emplace_back(fmt::format("{}{}{}", a, c, b));
    }
  }
}
}  // namespace xoos::demux
