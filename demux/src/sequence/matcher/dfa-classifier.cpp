#include "sequence/matcher/dfa-classifier.h"

#include <algorithm>
#include <queue>   // IWYU pragma: keep
#include <ranges>  // IWYU pragma: keep

#include "xoos/log/logging.h"
#include "xoos/util/sequence-functions.h"
#include "xxhash.h"

namespace xoos::demux {

namespace {

/**
 * @struct VectorHash
 * @brief Custom hasher for std::vector<s32> using xxHash during DFA compilation.
 */
struct VectorHash {
  std::size_t operator()(const std::vector<s32>& v) const noexcept {
    return XXH3_64bits(v.data(), v.size() * sizeof(s32));
  }
};

/// @brief Branchless C-locale ASCII uppercase. Avoids <cctype> functions (which depend on the global
/// locale) and the explicit-signedness char-to-int conversion they require.
constexpr char AsciiToUpper(const char c) {
  // Intentional ASCII character arithmetic (subtract to clear the lowercase bit) stored back as char.
  return (c >= 'a' && c <= 'z') ? static_cast<char>(c - ('a' - 'A')) : c;  // NOSONAR(S810)
}

}  // namespace

DFAClassifier::DFAClassifier(const std::string_view reference,
                             const SearchDirection direction,
                             const s32 max_drop_difference,
                             const s32 max_band_size)
    : _direction{direction} {
  // max_drop_difference can't be negative because we normalize the DFA to compute the drop and we'd prune every state.
  if (max_drop_difference < 0) {
    throw error::Error("DFAClassifier: max_drop_difference must be non-negative, got {}", max_drop_difference);
  }
  // The DFA needs a string larger than 0 in length
  if (reference.empty()) {
    throw error::Error("String for DFA construction was empty.");
  }

  // When direction is kBackward, reverse the reference so that Classify can iterate
  // the input sequence in reverse and still match against the correct pattern.
  const std::string reversed_storage(reference.rbegin(), reference.rend());
  const std::string_view effective_ref = direction == SearchDirection::kBackward ? reversed_storage : reference;

  // Map to track unique superpositions (NFA state sets) and assign them DFA index IDs
  // In NFA-to-DFA Powerset Construction, a single DFA state represents a "Set" of simultaneously active NFA states.
  // Here, the std::vector key is that mathematical set. The vector's size is exactly (effective_ref.size() + 1)
  // The index represents the position in the reference and the value is the log-odds score at that position.
  // For example:
  // A s32 vector of [0, -2, kSentinelPrunedScore, kSentinelPrunedScore] means we are most likely at Position 0,
  // possibly at Position 1 (with a -2 penalty), and paths to Position 2 and 3 are pruned.
  std::unordered_map<std::vector<s32>, u16, VectorHash> nfa_subset_to_dfa_idx;
  // keep track of NFA state sets
  std::queue<std::vector<s32>> unprocessed_nfa_subsets;

  // Initialize the start state superposition.
  auto start_scores = BuildStartSuperposition(effective_ref.size(), max_drop_difference);

  // state map key is a node/set (start_scores) with node index initialized to 0
  nfa_subset_to_dfa_idx[start_scores] = 0;
  unprocessed_nfa_subsets.emplace(start_scores);

  // Allocate State 0 in our vector (i.e. empty object)
  _states.emplace_back();
  // Start state: peak is at reference position 0, so all ref bases remain.
  _remaining_ref_bases.push_back(static_cast<u8>(effective_ref.size()));

  // Banded alignment depth limit: the maximum number of input bases the DFA can consume.
  // In banded alignment, the band constrains |input_pos - ref_pos| <= max_band_size.
  // A perfect alignment consumes ref_len input bases. With max_band_size extra insertions allowed,
  // the maximum input length is ref_len + max_band_size. States beyond this depth get pruned.
  const s32 max_depth = static_cast<s32>(effective_ref.size()) + max_band_size;

  // Track BFS depth for each state (depth = minimum number of input bases consumed to reach it)
  std::unordered_map<u16, s32> state_depth;
  state_depth[0] = 0;

  // Powerset Construction Loop (NFA -> DFA)
  // The NFA doesn't actually exist as an object but our superposition sets do
  while (!unprocessed_nfa_subsets.empty()) {
    auto current_subset = unprocessed_nfa_subsets.front();
    unprocessed_nfa_subsets.pop();

    const u16 current_dfa_idx = nfa_subset_to_dfa_idx[current_subset];
    const s32 current_depth = state_depth[current_dfa_idx];

    // Evaluate all 4 possible biological bases
    for (u8 base_idx = 0; base_idx < 4; ++base_idx) {
      const char base = sequence::IndexToBase(base_idx);
      // Advance the NFA superposition by one input base.
      // This is where the dynamic programming happens to simulate all possible paths.
      auto [next_scores, max_score] = AdvanceNfa(current_subset, base, effective_ref, max_drop_difference);

      // If all mathematical paths died, wire this transition to the pruned tombstone.
      if (max_score == kSentinelPrunedScore) {
        _states[current_dfa_idx].transitions[base_idx].PackSpecial(Tombstone::kPruneTermination, 0);
        continue;
      }

      // If the dominant path (score == 0) has reached the very end of the reference (L),
      // we have found the best possible biological match. Wire to the Success tombstone.
      if (next_scores[effective_ref.size()] == 0) {
        // Note: the cast to s8 is safe because max_score is mathematically bounded to [-4, 2]
        _states[current_dfa_idx].transitions[base_idx].PackSpecial(Tombstone::kSequenceEndMatch,
                                                                   static_cast<s8>(max_score));
        continue;
      }

      // Banded alignment depth check: if advancing would exceed the maximum allowed input bases,
      // wire to kPruneTermination. The alignment has consumed too many input bases without
      // completing the reference match, meaning it has exceeded the band.
      if (current_depth + 1 > max_depth) {
        _states[current_dfa_idx].transitions[base_idx].PackSpecial(Tombstone::kPruneTermination,
                                                                   static_cast<s8>(max_score));
        continue;
      }

      // Look up or register the target DFA state for this superposition
      const auto it = nfa_subset_to_dfa_idx.find(next_scores);
      u16 target_dfa_idx = 0;

      if (it != nfa_subset_to_dfa_idx.end()) {
        target_dfa_idx = it->second;
      } else {
        target_dfa_idx = static_cast<u16>(_states.size());

        if (target_dfa_idx >= static_cast<u16>(Tombstone::kSequenceEndMatch)) {
          throw error::Error(
              "DFAClassifier: State space exploded past 13-bit limit for sequence:\n{}\nRemaining unprocessed "
              "nfa subsets: {}. Try a smaller reference or stricter max_drop_difference. Max index "
              "possible: {}.",
              reference,
              unprocessed_nfa_subsets.size(),
              static_cast<u16>(Tombstone::kSequenceEndMatch) - 1);
        }

        nfa_subset_to_dfa_idx[next_scores] = target_dfa_idx;
        _states.emplace_back();

        // Compute remaining reference bases for this state.
        // The superposition is normalized so peak positions have score 0.
        // Find the furthest reference position at the peak to determine how
        // many reference bases the best NFA path still needs to cover.
        u8 furthest_peak = 0;
        for (size_t i = 0; i <= effective_ref.size(); ++i) {
          if (next_scores[i] == 0) {
            furthest_peak = static_cast<u8>(i);
          }
        }
        _remaining_ref_bases.push_back(static_cast<u8>(effective_ref.size() - furthest_peak));

        unprocessed_nfa_subsets.push(next_scores);
        state_depth[target_dfa_idx] = current_depth + 1;
      }

      // Wire the current state to the target, baking in the delta
      _states[current_dfa_idx].transitions[base_idx].Pack(target_dfa_idx, static_cast<s8>(max_score));
    }
  }
  // Log out some debugging
  Logging::Debug("DFA for {} created. {} direction. Max drop difference: {} Max band size: {}",
                 reference,
                 _direction == SearchDirection::kForward ? "Forward" : "Reverse",
                 max_drop_difference,
                 max_band_size);
  LogDebugStats();
}

std::vector<s32> DFAClassifier::BuildStartSuperposition(const size_t ref_len, const s32 max_drop_difference) {
  std::vector start_scores(ref_len + 1, kSentinelPrunedScore);
  start_scores[0] = 0;

  // Allow deletions right off the bat (sweeping left-to-right)
  for (size_t i = 0; i < ref_len; ++i) {
    const auto score = start_scores[i] + scoring::kDeletion;
    start_scores[i + 1] = score < -max_drop_difference ? kSentinelPrunedScore : score;
  }

  return start_scores;
}

DFAClassifier::NfaStepResult DFAClassifier::AdvanceNfa(const std::vector<s32>& current_subset,
                                                       const char base,
                                                       const std::string_view effective_ref,
                                                       const s32 max_drop_difference) {
  std::vector next_scores(effective_ref.size() + 1, kSentinelPrunedScore);

  // Simulates the NFA via dynamic programming
  for (size_t i = 0; i <= effective_ref.size(); ++i) {
    if (current_subset[i] == kSentinelPrunedScore) {
      continue;
    }

    // Insertion (The read advances, but the reference position stays the same)
    next_scores[i] = std::max(next_scores[i], current_subset[i] + scoring::kInsertion);

    // Matches & Substitutions (Both read and reference advance)
    if (i < effective_ref.size()) {
      const auto step_score = (base == AsciiToUpper(effective_ref[i])) ? scoring::kMatch : scoring::kSubstitution;
      next_scores[i + 1] = std::max(next_scores[i + 1], current_subset[i] + step_score);
    }
  }

  // Deletions (The reference advances, but the read stays the same)
  // This requires a left-to-right sweep so multiple consecutive deletions accumulate correctly.
  for (size_t i = 0; i < effective_ref.size(); ++i) {
    if (next_scores[i] != kSentinelPrunedScore) {
      next_scores[i + 1] = std::max(next_scores[i + 1], next_scores[i] + scoring::kDeletion);
    }
  }

  // Find the peak score across all reference positions
  auto max_score = kSentinelPrunedScore;
  for (const auto score : next_scores) {
    max_score = std::max(score, max_score);
  }

  // Normalize the vector so the highest score is always 0.
  // This is the mathematical trick that prevents infinite state explosion.
  if (max_score != kSentinelPrunedScore) {
    for (auto& score : next_scores) {
      if (score != kSentinelPrunedScore) {
        score -= max_score;
        // Apply X-drop pruning: prune paths that have dropped more than max_drop_difference from the peak
        if (score < -max_drop_difference) {
          score = kSentinelPrunedScore;
        }
      }
    }
  }

  return {std::move(next_scores), max_score};
}

void DFAClassifier::LogDebugStats() const {
  u32 pruned_count = 0;
  u32 match_count = 0;
  u32 self_loop_count = 0;
  const u32 total_transitions = static_cast<u32>(_states.size()) * 4;

  for (size_t s = 0; s < _states.size(); ++s) {
    for (const auto& t : _states[s].transitions) {
      const auto next = t.GetNextState();
      if (next == static_cast<u16>(Tombstone::kPruneTermination)) {
        ++pruned_count;
      } else if (next == static_cast<u16>(Tombstone::kSequenceEndMatch)) {
        ++match_count;
      } else if (next == static_cast<u16>(s)) {
        ++self_loop_count;
      }
    }
  }

  // BFS to find the maximum depth (longest shortest-path from state 0 to any tombstone)
  s32 max_depth = 0;
  {
    std::vector depth(_states.size(), -1);
    std::queue<u16> bfs_queue;
    depth[0] = 0;
    bfs_queue.push(0);
    while (!bfs_queue.empty()) {
      const auto state = bfs_queue.front();
      bfs_queue.pop();
      for (const auto& t : _states[state].transitions) {
        if (const auto next = t.GetNextState(); next < _states.size() && depth[next] == -1) {
          depth[next] = depth[state] + 1;
          max_depth = std::max(max_depth, depth[next]);
          bfs_queue.push(next);
        }
      }
    }
  }

  Logging::Debug("DFA has {} states, {} transitions ({} pruned, {} match, {} self-loop), max depth {}, {} bytes",
                 _states.size(),
                 total_transitions,
                 pruned_count,
                 match_count,
                 self_loop_count,
                 max_depth,
                 _states.size() * sizeof(DFAState));
}

template <typename Range>
DFAClassifier::DFAResults DFAClassifier::ScoreWalk(Range range, const s32 x_drop) const {
  using State = DFAResult::State;

  u16 current_state = _start_state;
  s32 current_score = 0;
  s32 max_seen_score = 0;

  DFAResult best{0, 0, State::kUnresolved};
  std::optional<DFAResult> last_match;
  // Track the DFA state we transitioned to on the last genuine base match (delta > 0).
  // Initialised to _start_state so that "no match ever seen" resolves to the full reference length when we look up
  // _remaining_ref_bases once after the loop.
  u16 last_match_next_state = _start_state;
  s32 pos = 0;
  bool exited_early = false;

  // The classification inner loop.
  // Designed for extreme cache-locality and zero complex branching.
  for (const char c : range) {
    // increment position. pos = 0 refers to no movement
    ++pos;
    const u8 base_idx = sequence::BaseToIndex(c);
    const Transition& t = _states[current_state].transitions[base_idx];

    // Always apply the delta immediately
    const s8 delta = t.GetScoreDelta();
    current_score += delta;

    // Resolve new state
    current_state = t.GetNextState();

    if (current_score > max_seen_score) {
      max_seen_score = current_score;
      best = {pos, current_score, State::kUnresolved};
    }

    // Track the last position where a genuine base match occurred (positive delta).
    // Only kMatch (+2) produces a positive delta; substitutions, insertions, and deletions are all negative.
    if (delta > 0) {
      last_match = DFAResult{pos, current_score, State::kUnresolved};
      last_match_next_state = current_state;
    }

    if (current_state >= static_cast<u16>(Tombstone::kSequenceEndMatch)) {
      const auto state =
          current_state == static_cast<u16>(Tombstone::kSequenceEndMatch) ? State::kFoundEnd : State::kFoundPruned;
      ResolveTerminalState(best, last_match, pos, state);
      exited_early = true;
      break;
    }

    // Runtime X-drop: terminate when cumulative score drops too far below the peak.
    // This is the BLAST-style X-drop — once the alignment has clearly degraded past the
    // adapter region, there is no value in continuing. All useful information (best,
    // last_match) has already been captured.
    if (max_seen_score - current_score > x_drop) {
      exited_early = true;
      break;
    }
  }

  // Reached end of input sequence without hitting a tombstone or X-drop
  if (!exited_early) {
    ResolveTerminalState(best, last_match, pos, State::kInputEnd);
  }

  return {.best = best,
          .last_match = last_match,
          .remaining_last_match_bases = ResolveRemainingBases(last_match_next_state)};
}

void DFAClassifier::ResolveTerminalState(DFAResult& best,
                                         std::optional<DFAResult>& last_match,
                                         const s32 pos,
                                         const DFAResult::State state) {
  if (best.pos == pos) {
    best.state = state;
  }
  // NOSONAR(M23_002) -- last_match->pos equals pos only when a positive-delta match occurred on this
  // iteration; the analyzer cannot see that pos advances independently, so both branches are reachable.
  if (last_match.has_value() && last_match->pos == pos) {  // NOSONAR(M23_002)
    last_match->state = state;
  }
}

s32 DFAClassifier::ResolveRemainingBases(const u16 last_match_next_state) const {
  // Single access to _remaining_ref_bases after the hot path (i.e. L1 cache contains DFA).
  // If last_match_next_state is a tombstone (only possible when the last match itself triggered a
  // tombstone on the same step), resolve without a table lookup:
  //   kSequenceEndMatch -> 0 (reference fully consumed)
  //   kPruneTermination -> fall back to start state (conservative; only reachable via the banded
  //                        depth limit and very rare)
  if (last_match_next_state == static_cast<u16>(Tombstone::kSequenceEndMatch)) {
    return 0;
  }
  if (last_match_next_state > static_cast<u16>(Tombstone::kSequenceEndMatch)) {
    return static_cast<s32>(_remaining_ref_bases[_start_state]);
  }
  return static_cast<s32>(_remaining_ref_bases[last_match_next_state]);
}

DFAClassifier::DFAResults DFAClassifier::Classify(const std::string_view sequence, const s32 x_drop) const {
  return _direction == SearchDirection::kBackward ? ScoreWalk(std::ranges::reverse_view(sequence), x_drop)
                                                  : ScoreWalk(sequence, x_drop);
}

}  // namespace xoos::demux
