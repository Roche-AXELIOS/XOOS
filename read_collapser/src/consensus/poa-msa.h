#pragma once

#include <string>

#include <spoa/spoa.hpp>

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::read_collapser {

/**
 * A multiple sequence aligner using the Partial Order Alignment (POA) algorithm.
 */
class PoaAligner {
 public:
  /**
   * Constructor to initialize the POA aligner with scoring parameters.
   *
   * @param match_score The score for a match. Must be positive.
   * @param mismatch_penalty The penalty for a mismatch. Must be non-negative.
   * @param gap_open_penalty The penalty for opening a gap. Must be non-negative.
   * @param gap_ext_penalty The penalty for extending a gap. Must be non-negative.
   */
  explicit PoaAligner(s32 match_score, s32 mismatch_penalty, s32 gap_open_penalty, s32 gap_ext_penalty);

  /**
   * Given a list of sequences, compute the multiple sequence alignment (MSA)
   * using the partial order alignment (POA) algorithm. Partial order alignment
   * creates a multiple sequence alignment by progressively aligning sequences
   * onto a directed acyclic graph where each node represents a base and edges
   * represent adjacency between bases. Once the graph is built, the MSA can be
   * extracted by traversing the graph in topological order.
   *
   * @param sequences The input sequences to align.
   * @return A vector of aligned sequences.
   */
  vec<std::string> ComputeMSA(const vec<std::string>& sequences);

 private:
  std::unique_ptr<spoa::AlignmentEngine> _engine;
};

}  // namespace xoos::read_collapser
