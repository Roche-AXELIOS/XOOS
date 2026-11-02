#pragma once

#include <optional>
#include <string_view>
#include <vector>

#include <xoos/error/error.h>
#include <xoos/types/int.h>

#include "sequence/alignment/log-likelihood-scoring.h"
#include "sequence/matcher/search-direction.h"

namespace xoos::demux {

/**
 * @class DFAClassifier
 * @brief  This is a Deterministic Finite Automaton (DFA)-based sequence classifier, which is designed to efficiently
 * classify sequences based on a single reference pattern (e.g., adapter stem), for seed-and-extend approaches.
 *
 * This is originally intended when we know the start of where the match should be. Matching ends when we reach the end
 * of the reference sequence not penalizing for trailing insertions or even missing sequence entirely with k-mer based
 * methods (because of truncated sequence).
 *
 * Unlike purely Levenshtein-based approaches, such as Myers' Bit-Vector Algorithm, this DFA-based classifier can
 * incorporate more complex heuristics and scoring schemes, such as different penalties for insertions, deletions, and
 * substitutions. This allows for a more nuanced classification that can better reflect the biological realities of
 * sequencing errors and variations. The downside is that the DFA can grow exponentially with the length of the
 * reference sequence, so it is important to optimize the construction and prune the state space effectively.
 *
 * Technically we are performing "Powerset Construction" of an NFA to a DFA, but we can optimize by bypassing the
 * explicit NFA construction. We do this by directly building the DFA states and transitions based on the reference
 * sequence and the defined scoring scheme for matches, insertions, substitutions, and deletions
 *
 * The classifier uses log odds values for scoring matches, insertions, substitutions, and deletions. These values are
 * additive for fast computation and can be tuned based on empirical data to reflect the likelihood of different types
 * of errors in the sequencing process. Analogous to scoring and will not be a true probability, but a log odds score
 * that can be used for classification thresholds. Could be calibrated to actual probabilities if we ever want that.
 *
 * Pruning during construction:
 * - NFA superposition width limit (max_drop_difference): during powerset construction, each NFA superposition
 *   vector is normalized so the peak position has score 0. Individual NFA positions whose normalized score falls
 *   below -max_drop_difference are pruned. This controls the number of distinct DFA states (memory) and limits
 *   the error patterns the DFA can represent. For example, max_drop_difference=10 allows up to 2 consecutive
 *   deletions (2 * |kDeletion| = 10) within a single alignment step to survive.
 *   Note: Due to the insertion self-loops (the peak position can always absorb an input base via insertion at
 *   cost kInsertion, which after normalization can never exceed -|kMatch - kInsertion| = -5), so we technically can't
 *   ever end on these nodes in a fully pruned state in the graph.
 *
 * - Banded alignment depth limit (max_band_size): During powerset construction, the BFS depth of each state
 *   is tracked (depth = minimum number of input bases consumed to reach it). The maximum allowed depth is
 *   ref_len + max_band_size, analogous to banded Smith-Waterman where the alignment can deviate at most
 *   max_band_size positions from the diagonal. Any transition that would create a state beyond this depth
 *   is wired to kPruneTermination instead. This bounds the DFA size and produces structural dead-ends
 *   that terminate classification when too many input bases are consumed without completing the reference.
 *
 * Pruning during classification (runtime):
 * - Runtime X-drop in Classify: when the cumulative score drops more than a caller-specified threshold below the
 *   peak score observed so far, the walk terminates early with kFoundPruned. This is the BLAST-style X-drop and
 *   is the primary mechanism for stopping once the alignment has clearly moved past the adapter region.
 */
class DFAClassifier {
 public:
  struct DFAResult {
    enum class State {
      // state where we could still continue but ended for some reason like an X-drop threshold
      kUnresolved,
      // Full reference was consumed — clean alignment endpoint
      kFoundEnd,
      // Hit a structural kPruneTermination tombstone during the DFA walk
      kFoundPruned,
      // Exhausted the input sequence without hitting a tombstone
      kInputEnd,
    };

    // Position within the input sequence where the score was observed
    s32 pos;
    // Score observed during the walk by this state
    s32 score;
    // the state type we saw at this position
    State state;
  };

  struct DFAResults {
    // Best state we ended on, useful for allowing for flank after sequence not to penalize match
    DFAResult best{};
    // Last match state we ended on, optional because it is possible we saw no matches
    // Useful for when we classification the match to read teh end of seqeunces (e.g. trim boundries)
    std::optional<DFAResult> last_match;
    // Minimum number of reference positions still unmatched on the best NFA path at the point the walk terminated.
    // Can be up to the sequence length if we didn't have a last match
    s32 remaining_last_match_bases{};
  };

  /**
   * Constructs the DFA classifier with the given reference pattern.
   *
   * We construct the DFA by directly building the states and transitions based on the reference sequence.
   * When direction is kBackward the reference is reversed internally so that Classify iterates the input
   * sequence in reverse, producing biologically meaningful results without the caller needing to reverse anything.
   *
   * The max_drop_difference controls the width of the NFA superposition during powerset construction. The NFA
   * superposition is a vector of scores representing the best log-odds score to reach each position in the reference at
   * a given point in the input sequence. After each AdvanceNfa step, the superposition is normalized so the peak
   * position has score 0. Any NFA position whose normalized score falls below -max_drop_difference is pruned from the
   * superposition.
   *
   * Banded alignment width determines how far out indels can go out to. For example, to tolerate 2 consecutive
   * deletions in an alignment, this must be >= 2 * |kDeletion| = 10. The DFA will allow at most ref_len + max_band_size
   * input bases before pruning.
   *
   * @param reference The reference sequence to build the automaton against
   * @param direction The intended scanning direction (kForward or kBackward). Defaults to kForward.
   * @param max_drop_difference Controls the width of the NFA superposition during powerset construction.
   * @param max_band_size Banded alignment width, analogous to the bandwidth in banded Smith-Waterman.
   */
  explicit DFAClassifier(std::string_view reference,
                         SearchDirection direction = SearchDirection::kForward,
                         s32 max_drop_difference = 5,
                         s32 max_band_size = 10);

  /**
   * Classifies a given sequence against the compiled automaton in the direction specified at construction time.
   *
   * For kForward: scans from the beginning of the sequence towards the end.
   * For kBackward: scans from the end of the sequence towards the beginning.
   *
   * Terminates when:
   * - The DFA reaches a tombstone state (kSequenceEndMatch or kPruneTermination) -> kFoundPruned
   * - The cumulative score drops more than x_drop below the peak score (runtime X-drop) -> kUnresolved
   * - The input sequence is exhausted -> kInputEnd
   *
   * @param sequence The sequence to classify
   * @param x_drop  Runtime X-drop threshold. If the cumulative score falls more than x_drop below the highest score
   *                observed so far, the walk terminates early.
   * @return DFAResults containing best and last_match results from the walk.
   */
  DFAResults Classify(std::string_view sequence, s32 x_drop = 12) const;

 private:
  /**
   * Shared inner walk for Classify. The forward and backward paths differ only in the iteration order
   * over the input, so both instantiate this template with the appropriate range view. Kept as a named
   * member (rather than a lambda inside Classify) so its control flow is analyzed on its own.
   *
   * @tparam Range A cheap range view of char yielding the input in the desired scan order.
   * @param range The input sequence in scan order (a lightweight view, taken by value).
   * @param x_drop Runtime X-drop threshold (see Classify).
   * @return DFAResults for this walk.
   */
  template <typename Range>
  DFAResults ScoreWalk(Range range, s32 x_drop) const;

  /**
   * Stamp the terminal State onto best and/or last_match when they end at the current position.
   * Shared by the tombstone-exit and input-exhausted paths of ScoreWalk.
   */
  static void ResolveTerminalState(DFAResult& best,
                                   std::optional<DFAResult>& last_match,
                                   s32 pos,
                                   DFAResult::State state);

  /**
   * Resolve the number of reference bases still outstanding after the last genuine match, given the
   * DFA state reached on that match. Tombstone states are handled without a table lookup.
   */
  s32 ResolveRemainingBases(u16 last_match_next_state) const;

  /// Result of advancing the NFA superposition by one input base.
  struct NfaStepResult {
    std::vector<s32> next_scores;
    /// The unnormalized peak score across all reference positions before normalization.
    /// If all paths died this will be kSentinelPrunedScore and next_scores is meaningless.
    s32 max_score;
  };

  /**
   * Advance the NFA superposition by one input base.
   *
   * Given the current superposition (score at every reference position) and an input base,
   * computes the next superposition by simulating insertion, match/substitution, and deletion
   * transitions via dynamic programming, then normalizes and X-drop prunes the result.
   *
   * @param current_subset  Current NFA superposition vector of size (ref_len + 1)
   * @param base            The input base character (A, C, G, or T)
   * @param effective_ref   The reference sequence being matched against
   * @param max_drop_difference Pruning from normalize peak threshold
   * @return NfaStepResult containing the normalized next superposition and the unnormalized peak score
   */
  static NfaStepResult AdvanceNfa(const std::vector<s32>& current_subset,
                                  char base,
                                  std::string_view effective_ref,
                                  s32 max_drop_difference);

  /**
   * Builds the initial NFA superposition vector for the start state.
   *
   * Position 0 gets score 0; subsequent positions are filled by sweeping deletions left-to-right,
   * pruning any that fall below -max_drop_difference.
   *
   * @param ref_len            Length of the effective reference sequence
   * @param max_drop_difference  Pruning threshold (positive value)
   * @return Superposition vector of size (ref_len + 1)
   */
  static std::vector<s32> BuildStartSuperposition(size_t ref_len, s32 max_drop_difference);

  /**
   * Logs debug statistics about the DFA, such as the number of states, transitions, and distribution of score deltas.
   */
  void LogDebugStats() const;

  // Sentinel value representing a dead/pruned mathematical path in the NFA superposition.
  // Chosen to be far from any real score (which are normalized to [−max_drop_difference, 0]) but not so extreme
  // as to cause signed overflow when penalties are added to it.
  static constexpr s32 kSentinelPrunedScore = -99999;

  // Bitwise operations on signed negative numbers can cause sign-extension bugs
  // Our delta range is [kDeletionPenalty, kMatchScore] = [-5, 2]. Adding offset of +5 maps to [0, 7].
  static constexpr s8 kDeltaOffset = -scoring::kMinPenalty;
  // 3 bits (7 values needed total from -4 to 2) for the score delta
  static constexpr auto kDeltaBits = 3;
  // Mask to extract the score delta (3 bits)
  static constexpr u16 kDeltaMask = (1 << kDeltaBits) - 1;
  // Remaining 13 bits (2^13 - 1 = 8191) are used to store indexes that we can jump to
  // Mask to extract the index (13 bits)
  static constexpr u16 kIndexMask = (1 << (16 - kDeltaBits)) - 1;

  // Tombstone flags encoded in the 13-bit index field, replacing normal state indexes with special meanings.
  // Values at the top of the index range are reserved as tombstones for control flow during classification.
  enum class Tombstone : u16 {
    // This tombstone (8191) indicates that we should terminate the search because it has been pruned
    kPruneTermination = (1 << (16 - kDeltaBits)) - 1,
    // This tombstone (8190) indicates that we have reached an ending that satisfies the reference sequence
    kSequenceEndMatch = (1 << (16 - kDeltaBits)) - 2
  };

  /**
   * @struct Transition
   * @brief Represents a transition in the DFA, encoding both the next state index and the score delta.
   *
   * The `packed_data` member is a 16-bit unsigned integer that encodes both the next state index (13 bits) and the
   * score delta (3 bits). The `Pack` method is used to pack these two pieces of information into a single integer,
   * while the `GetNextState` and `GetScoreDelta` methods are used to extract the next state index and score delta
   * from the packed data, respectively.
   */
  struct Transition {
    u16 packed_data;

    /**
     * Packs the next state index and score delta into a single 16-bit unsigned integer.
     *
     * @param next_state_index The index of the next state to transition to (must be less than 8190)
     * @param score_delta The change in score associated with this transition (must be between -4 and +2)
     * @throws error::Error if next_state_index is out of bounds or if score_delta is out of range
     */
    void Pack(u16 next_state_index, const s8 score_delta) {
      if (next_state_index >= static_cast<u16>(Tombstone::kSequenceEndMatch)) {
        throw error::Error("DFAClassifier::Transition::Pack: next_state_index must be less than {}. Got {}.",
                           static_cast<u16>(Tombstone::kSequenceEndMatch),
                           next_state_index);
      }
      // Shift state left by 3 bits, then insert the offset delta into the bottom 3 bits. Compute in
      // u32 so the u16 operands are not integer-promoted to (signed) int (M23_380).
      const auto unsigned_delta = static_cast<u16>(score_delta + kDeltaOffset);
      const u32 index_bits = static_cast<u32>(next_state_index) << static_cast<u32>(kDeltaBits);
      const u32 delta_bits = static_cast<u32>(unsigned_delta) & static_cast<u32>(kDeltaMask);
      packed_data = static_cast<u16>(index_bits | delta_bits);
    }

    /**
     * Packs a special state flag into the transition, which indicates a special control flow (e.g., termination).
     *
     * @param special_flag The tombstone flag to pack (must be a valid Tombstone)
     * @param score_delta The change in score associated with this transition (must be between -4 and +2)
     * @throws error::Error if special_flag is not a valid Tombstone or if score_delta is out of range
     */
    void PackSpecial(Tombstone special_flag, const s8 score_delta) {
      using enum Tombstone;
      switch (special_flag) {
        case kPruneTermination:
        case kSequenceEndMatch:
          // continue as normal
          break;
        default:
          throw error::Error("DFAClassifier::Transition::PackSpecial: special_flag must be a valid Tombstone. Got {}.",
                             static_cast<u16>(special_flag));
      }
      // Compute in u32 so the u16 operands are not integer-promoted to (signed) int (M23_380).
      const auto unsigned_delta = static_cast<u16>(score_delta + kDeltaOffset);
      const u32 flag_bits = static_cast<u32>(special_flag) << static_cast<u32>(kDeltaBits);
      const u32 delta_bits = static_cast<u32>(unsigned_delta) & static_cast<u32>(kDeltaMask);
      packed_data = static_cast<u16>(flag_bits | delta_bits);
    }

    // Inline helpers for the runtime loop
    u16 GetNextState() const {
      // Shift right to drop the delta bits
      return packed_data >> kDeltaBits;
    }

    s8 GetScoreDelta() const {
      return static_cast<s8>((packed_data & kDeltaMask) - kDeltaOffset);
    }
  };

  // 16 * 4 bits = 64 bits per state
  struct DFAState {
    Transition transitions[4];
  };

  SearchDirection _direction;
  u16 _start_state = 0;
  // The actual states used
  std::vector<DFAState> _states;
  // Parallel to _states: for each DFA state, the minimum number of reference bases remaining on the best NFA path
  std::vector<u8> _remaining_ref_bases;
};

}  // namespace xoos::demux
