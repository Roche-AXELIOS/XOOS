#include "consensus/poa-msa.h"

#include <spoa/spoa.hpp>

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "consensus/base-encoder.h"

namespace xoos::read_collapser {

PoaAligner::PoaAligner(const s32 match_score,
                       const s32 mismatch_penalty,
                       const s32 gap_open_penalty,
                       const s32 gap_ext_penalty) {
  _engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW,
                                          static_cast<s8>(match_score),
                                          static_cast<s8>(-mismatch_penalty),
                                          static_cast<s8>(-gap_open_penalty),
                                          static_cast<s8>(-gap_ext_penalty));
}

vec<std::string> PoaAligner::ComputeMSA(const vec<std::string>& sequences) {
  // Create a SPOA graph
  spoa::Graph graph;

  // Add sequences to the graph
  for (const auto& sequence : sequences) {
    if (sequence.empty()) {
      // SPOA does not align empty sequences, so we need to handle them separately.
      continue;
    } else {
      // For non-empty sequences, align and add to the graph
      graph.AddAlignment(_engine->Align(sequence, graph), sequence);
    }
  }

  // Convert the alignment to a vector of strings
  auto msa = graph.GenerateMultipleSequenceAlignment();

  if (msa.empty()) {
    // All sequences are empty so we return a vector of gaps
    return {sequences.size(), std::string(1, kBaseGap)};
  } else {
    // Some sequences are non-empty, we need to merge the MSA result with empty sequences
    vec<std::string> msa_with_empty_sequences;
    size_t j = 0;
    for (const auto& sequence : sequences) {
      if (sequence.empty()) {
        // Since SPOA does not align empty sequences, they are not part of the MSA result.
        // We need to fill them in manually with gap sequences.
        // Since at this point we know that the msa is not empty and that all sequences in the
        // msa result have the same length, we can use the length of the first aligned sequence
        // as the length of the gap sequence.
        msa_with_empty_sequences.emplace_back(msa[0].length(), kBaseGap);
      } else {
        // We won't go out of bounds because msa.size() <= sequences.size()
        // and we only increment j when we encounter a non-empty sequence
        msa_with_empty_sequences.emplace_back(msa[j]);
        ++j;
      }
    }
    return msa_with_empty_sequences;
  }
}

}  // namespace xoos::read_collapser
