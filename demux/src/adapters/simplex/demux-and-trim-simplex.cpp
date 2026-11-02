#include "demux-and-trim-simplex.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <string_view>
#include <utility>

#include <magic_enum/magic_enum.hpp>

#include "xoos/log/logging.h"
#include "xoos/types/vec.h"

namespace xoos::demux {

namespace {
// Upper bound on the anchor length fed to an extension DFA.
//
// The DSL config (source of truth) can carry long fixed anchors (e.g. the 21bp trailing anchor
// TTTGTCGTGTAGGTGAAGAGC, or the 29-34bp anchors of the I2X-S-DX dual-index design). An extension DFA
// only bounds how far trimming reaches from the SID boundary into the anchor, so it does not need the
// full anchor -- and it cannot take it: DFA powerset construction explodes past its 13-bit state
// limit once the reference is long enough at the configured max_drop_difference. The explosion length
// is content- and direction-dependent, so no single constant is correct for every design.
//
// This bound (15) is chosen to be a no-op for every current production simplex design (YS anchors are
// <= 15bp, with fixed_7 historically pinned to a 15bp prefix), preserving their behaviour exactly.
// MakeExtendDfa then backs off further, one base at a time, only when construction still explodes
// for a given anchor/direction/drop combination (as the longer I2X anchors require).
constexpr std::size_t kExtendAnchorDfaMaxLen = 15;

// Upper bound on the stem length fed to an extension DFA.
//
// Unlike the outer anchors, the stems sit between the SID and the insert, so the extension DFA must
// span the WHOLE stem to trim it off the insert -- any bases it cannot reach are left attached to the
// demuxed read. There is therefore no legitimate cap: this bound is set high enough to never truncate
// a real stem (the longest current design, I2X-S-DX, has 33-34bp stems). MakeExtendDfa still backs
// off per base if a stem genuinely explodes the DFA state space, but at the production max_drop the
// full stems build cleanly, so backoff never triggers for current designs.
constexpr std::size_t kExtendStemDfaMaxLen = 1024;

// Read-order positions of the simplex adapter components, named locally at the single point where
// the components are addressed. The layout is [anchor, sid, stem, stem, sid, anchor] flanking the
// insert: 5' [anchor_fwd][sid_fwd][stem_fwd] insert [stem_rev][sid_rev][anchor_rev] 3'. These index
// a component vector directly, so they are plain constants rather than a scoped enum.
constexpr std::size_t kAnchorFwd = 0U;
constexpr std::size_t kSidFwd = 1U;
constexpr std::size_t kStemFwd = 2U;
constexpr std::size_t kStemRev = 3U;
constexpr std::size_t kSidRev = 4U;
constexpr std::size_t kAnchorRev = 5U;

// max_drop_difference for each extension DFA: how far the alignment score may fall below its running
// peak before pruning during powerset construction (see DFAClassifier). Larger values tolerate more
// consecutive indels in the fixed region at the cost of a wider DFA. These are per-region and must
// stay in sync with the member-initializer order below.
//
// The stems use a smaller max_drop (4) than the anchors on purpose. The stems must be built at FULL
// length so the whole stem is trimmed off the insert (see kExtendStemDfaMaxLen and the
// stem-trimming bug note there). For the long I2X-S-DX stems (34bp / 33bp) the DFA powerset state
// count grows with both reference length and max_drop, and at max_drop >= 5 a full-length stem blows
// the 13-bit (8191) state limit in DFAClassifier -- MakeExtendDfa would then silently back off to a
// ~10bp slice, leaving ~2/3 of the stem attached to the demuxed read. max_drop = 4 is the largest
// value at which both full-length stems build cleanly (verified empirically), while still tolerating
// a couple of mismatches / one indel across the stem. The anchors face the read ends, so they can
// keep the wider max_drop (5) and are length-capped instead (see kExtendAnchorDfaMaxLen).
// anchor_fwd, 5' side, backward
constexpr s32 kBwExtendAnchorFwdMaxDrop = 5;
// Both independent SID stems use the narrow band required by the long I2X-S-DX design.
constexpr s32 kIndependentStemMaxDrop = 4;
// anchor_rev, 3' side, forward
constexpr s32 kFwExtendAnchorRevMaxDrop = 5;

// One-SID completion searches only the adapter-proximal ends of a read. The I2X-S-DX
// 5' and 3' stems are respectively 34bp and 33bp, with up to 37bp of SID plus outer adapter
// between a stem and the physical read end. 128bp leaves substantial slack for partial outer
// adapters while keeping the search space bounded and independent of read length.
constexpr s32 kPartialStemEndpointWindow = 128;
constexpr s32 kPartialStemMaxEditDistance = 4;
constexpr std::size_t kMinPartialStemLength = 16U;
constexpr u32 kPartialStemMaxErrorRatePercent = 5U;
constexpr u32 kPercentDenominator = 100U;

using PartialStemBitap = Bitap<kPartialStemMaxEditDistance>;

struct PartialStemCandidate {
  LociRange range;
  u32 edit_distance;
};

struct TerminalPartialStemCandidate {
  LociRange range;
  u32 expected_stem_length;
  u32 edit_distance;
};

struct BandedAlignment {
  std::string_view query;
  std::string_view observed;
};

struct BandedCell {
  std::size_t row;
  std::size_t column;
  std::size_t current_index;
};

struct BandedRow {
  const vec<u32>& previous;
  const vec<u32>& current;
  std::size_t previous_begin;
  std::size_t band_size;
  u32 infinity;

  bool ContainsPreviousColumn(const std::size_t column) const {
    return column >= previous_begin && column < previous_begin + band_size;
  }
};

u32 BoundedEditDistanceCell(const BandedAlignment& alignment, const BandedCell& cell, const BandedRow& row) {
  if (cell.row == 0U) {
    return static_cast<u32>(cell.column);
  }

  auto distance = row.infinity;
  if (cell.column > 0U && row.ContainsPreviousColumn(cell.column - 1U)) {
    distance = std::min(distance,
                        row.previous[cell.column - 1U - row.previous_begin] +
                            (alignment.query[cell.row - 1U] == alignment.observed[cell.column - 1U] ? 0U : 1U));
  }
  if (cell.current_index > 0U) {
    distance = std::min(distance, row.current[cell.current_index - 1U] + 1U);
  }
  if (row.ContainsPreviousColumn(cell.column)) {
    distance = std::min(distance, row.previous[cell.column - row.previous_begin] + 1U);
  }
  return std::min(distance, row.infinity);
}

// Validate the exact interval that will be trimmed. The two directional Bitap searches efficiently
// find candidate boundaries, but they are independent approximate searches and therefore cannot by
// themselves establish the edit distance of the interval formed by joining those boundaries.
std::optional<u32> BoundedEditDistance(const std::string_view query,
                                       const std::string_view observed,
                                       const u32 max_edit_distance) {
  const auto query_length = query.size();
  const auto observed_length = observed.size();
  const auto band_width = static_cast<std::size_t>(max_edit_distance);
  if (query_length > observed_length + band_width || observed_length > query_length + band_width) {
    return std::nullopt;
  }

  if (max_edit_distance == 0U) {
    return query == observed ? std::optional<u32>{0U} : std::nullopt;
  }

  const auto infinity = max_edit_distance + 1U;
  const auto band_size = (2U * band_width) + 1U;
  vec<u32> previous(band_size, infinity);
  vec<u32> current(band_size, infinity);
  const BandedAlignment alignment{.query = query, .observed = observed};

  // Only cells within max_edit_distance of the main diagonal can produce a result within the bound.
  // Rolling rows keep the memory use proportional to the band rather than to the full interval.
  for (std::size_t row = 0; row <= query_length; ++row) {
    const auto row_begin = row > band_width ? row - band_width : 0U;
    const auto row_end = std::min(observed_length, row + band_width);
    std::ranges::fill(current, infinity);
    const auto previous_begin = row > band_width + 1U ? row - band_width - 1U : 0U;
    const BandedRow previous_row{
        .previous = previous,
        .current = current,
        .previous_begin = previous_begin,
        .band_size = band_size,
        .infinity = infinity,
    };

    for (auto column = row_begin; column <= row_end; ++column) {
      const auto current_index = column - row_begin;
      const BandedCell cell{.row = row, .column = column, .current_index = current_index};
      current[current_index] = BoundedEditDistanceCell(alignment, cell, previous_row);
    }
    std::swap(previous, current);
  }

  const auto final_begin = query_length > band_width ? query_length - band_width : 0U;
  const auto final_index = observed_length - final_begin;
  const auto distance = previous[final_index];
  return distance <= max_edit_distance ? std::optional<u32>{distance} : std::nullopt;
}

u32 AllowedPartialStemEditDistance(const std::size_t expected_stem_length) {
  // Floor the rate rather than rounding up: a terminal fragment has no second adapter feature to
  // corroborate it, so the shorter lengths must be exact. This keeps accidental insert trimming
  // rare while allowing one error once at least 20 stem bases are observed.
  return static_cast<u32>((expected_stem_length * kPartialStemMaxErrorRatePercent) / kPercentDenominator);
}

std::optional<TerminalPartialStemCandidate> FindTerminalPartialStemCandidate(const std::string_view stem,
                                                                             const std::string_view read,
                                                                             const bool is_5p) {
  std::optional<TerminalPartialStemCandidate> best;

  // The acceptance of an approximate alignment is not monotonic with stem length, so binary
  // search is unsound. There are only 25 candidate lengths for current 33/34bp stems.
  for (std::size_t expected_length = stem.size(); expected_length >= kMinPartialStemLength; --expected_length) {
    const auto max_edit_distance = AllowedPartialStemEditDistance(expected_length);
    const auto min_observed_length =
        std::max(kMinPartialStemLength,
                 expected_length > max_edit_distance ? expected_length - max_edit_distance : std::size_t{0U});
    const auto max_observed_length = std::min(read.size(), expected_length + max_edit_distance);
    if (min_observed_length > max_observed_length) {
      continue;
    }

    const auto expected = is_5p ? stem.substr(stem.size() - expected_length) : stem.substr(0U, expected_length);
    for (std::size_t observed_length = min_observed_length; observed_length <= max_observed_length; ++observed_length) {
      const auto observed =
          is_5p ? read.substr(0U, observed_length) : read.substr(read.size() - observed_length, observed_length);
      const auto edit_distance = BoundedEditDistance(expected, observed, max_edit_distance);
      if (!edit_distance.has_value()) {
        continue;
      }

      const LociRange range =
          is_5p ? LociRange{0U, static_cast<u32>(observed_length)}
                : LociRange{static_cast<u32>(read.size() - observed_length), static_cast<u32>(read.size())};
      const TerminalPartialStemCandidate candidate{
          .range = range,
          .expected_stem_length = static_cast<u32>(expected_length),
          .edit_distance = *edit_distance,
      };
      if (!best.has_value() || candidate.expected_stem_length > best->expected_stem_length ||
          (candidate.expected_stem_length == best->expected_stem_length &&
           candidate.edit_distance < best->edit_distance)) {
        best = candidate;
      }
    }
  }
  return best;
}

bool HasPlausibleOuterContext(const std::string_view anchor,
                              const std::string_view outer_context,
                              const std::size_t unresolved_sid_length,
                              const bool is_5p) {
  if (outer_context.size() <= unresolved_sid_length) {
    return true;
  }

  return FindTerminalPartialStemCandidate(anchor, outer_context, is_5p).has_value();
}

std::optional<PartialStemCandidate> FindTerminalStemCandidate(const FixedReadRecord& record,
                                                              const bool is_5p,
                                                              const PartialStemBitap& forward,
                                                              const PartialStemBitap& backward,
                                                              const s32 window_begin,
                                                              const s32 window_end) {
  // A calibrated stem score can eventually replace these two directional Bitap boundary searches
  // to select both boundaries together and explicitly control false-positive trimming.
  const auto edge = (is_5p ? forward : backward).FindWidestMatch(record.Seq(), window_begin, window_end);
  if (!edge.Found()) {
    return std::nullopt;
  }

  const auto recovery_begin =
      is_5p ? std::max(window_begin, edge.position - forward.GetQueryLength() - kPartialStemMaxEditDistance + 1)
            : edge.position;
  const auto recovery_end =
      is_5p ? edge.position
            : std::min(window_end, edge.position + forward.GetQueryLength() + kPartialStemMaxEditDistance);
  const auto other_edge = (is_5p ? backward : forward).FindWidest(record.Seq(), recovery_begin, recovery_end);
  if (other_edge == kNoMatchPosition) {
    return std::nullopt;
  }

  const LociRange range = is_5p ? LociRange{static_cast<u32>(other_edge), static_cast<u32>(edge.position + 1)}
                                : LociRange{static_cast<u32>(edge.position), static_cast<u32>(other_edge + 1)};
  const auto max_edit_distance = AllowedPartialStemEditDistance(forward.Query().size());
  const auto interval_edit_distance =
      BoundedEditDistance(forward.Query(),
                          std::string_view(record.Seq() + range.spos, static_cast<std::size_t>(range.Length())),
                          max_edit_distance);
  if (!interval_edit_distance.has_value()) {
    return std::nullopt;
  }

  return PartialStemCandidate{.range = range, .edit_distance = *interval_edit_distance};
}

bool IsBetterPartialStemCandidate(const PartialStemCandidate& candidate,
                                  const PartialStemCandidate& best,
                                  const bool is_5p) {
  if (candidate.edit_distance != best.edit_distance) {
    return candidate.edit_distance < best.edit_distance;
  }
  return is_5p ? candidate.range.epos > best.range.epos : candidate.range.spos < best.range.spos;
}

// Build an extension DFA from the longest usable slice of a fixed region (anchor or stem).
//
// Trimming is measured relative to the SID match boundary, so the DFA must retain the bases adjacent
// to the SID: the prefix for a forward DFA, the suffix for a backward DFA. Starting from max_len,
// shorten the slice one base at a time (dropping from the SID-distal end) until DFA construction
// succeeds, so long fixed regions that would otherwise blow the powerset-state limit still yield a
// usable extension DFA instead of aborting the run.
//
// IMPORTANT: max_len must be chosen per region. For the OUTER anchors (anchor_fwd, anchor_rev) the
// SID-distal end faces the read end, so shortening (or capping via max_len) only limits how far
// trimming reaches into the anchor towards the read end -- never leaving adapter attached to the
// insert. Anchors therefore cap at kExtendAnchorDfaMaxLen to stay within the DFA state limit for the
// long I2X-S-DX anchors. For the INNER stems (stem_fwd, stem_rev) the SID-distal end faces the
// insert, so dropping those bases would leave the stem tail attached to the demuxed insert. Stems
// must therefore pass their full length as max_len (backoff still guards against a genuine state
// explosion, but at the production max_drop the full stems build cleanly).
DFAClassifier MakeExtendDfa(const std::string_view anchor,
                            const SearchDirection direction,
                            const s32 max_drop,
                            const std::size_t max_len) {
  std::size_t len = std::min(anchor.size(), max_len);
  while (len > 0) {
    const std::string_view slice =
        direction == SearchDirection::kBackward ? anchor.substr(anchor.size() - len) : anchor.substr(0, len);
    try {
      return DFAClassifier{slice, direction, max_drop};
    } catch (const std::runtime_error& e) {  // NOSONAR(S1181,S2486) DFAClassifier throws exactly std::runtime_error
      // (via error::Error), so no more-specific type exists; the state-space explosion is handled here
      // by logging and retrying with a shorter slice.
      Logging::Debug("Extension DFA for region '{}' exceeded the state limit at length {}; shortening. ({})",
                     anchor,
                     len,
                     e.what());
    }
    --len;
  }
  throw error::Error("Unable to build extension DFA for region '{}': state space explodes at every length", anchor);
}

// Running (score, trim position) for a candidate as an extension DFA's terminal state is folded in.
struct ScoreTrim {
  s32 score;
  s32 trim_pos;
};

// Fold a forward extension DFA's terminal state into a 5' candidate's score and trim position.
//
// The forward DFA walks the fixed region right of the SID, towards the insert, so any bases matched
// beyond the SID push the 5' trim boundary deeper (trim_pos increases). The terminal state decides how
// the tail is scored: a clean end adds its score as-is; an X-drop/unresolved end treats the remaining
// reference bases as substitutions and extends the trim over them; a pruned/input-end tail treats them
// as deletions (which do not consume read bases, so the trim boundary is left where it is).
ScoreTrim FoldForwardExtension(const DFAClassifier::DFAResults& fw, s32 score, s32 trim_pos) {
  using enum DFAClassifier::DFAResult::State;
  if (!fw.last_match.has_value()) {
    // No match anywhere in the forward region; charge every reference base as a substitution.
    const auto remaining_bases = fw.remaining_last_match_bases;
    return {score + (remaining_bases * scoring::kSubstitution), trim_pos + remaining_bases};
  }
  const auto& last = fw.last_match.value();
  // Add the bases matched beyond the SID to the trim boundary.
  trim_pos += last.pos;
  const auto remaining_bases = fw.remaining_last_match_bases;
  switch (last.state) {
    case kFoundEnd:
      score += last.score;
      break;
    case kUnresolved:
      score += last.score + (remaining_bases * scoring::kSubstitution);
      trim_pos += remaining_bases;
      break;
    case kFoundPruned:
    case kInputEnd:
      // Ran off the (likely truncated) end; the unmatched reference bases must be deletions.
      score += last.score + (remaining_bases * scoring::kDeletion);
      break;
    default:
      throw error::Error("Found impossible DFA state end {} ", magic_enum::enum_name(last.state));
  }
  return {score, trim_pos};
}

// Fold a backward extension DFA's terminal state into a 3' candidate's score and trim position.
//
// Mirror image of FoldForwardExtension: the backward DFA walks the fixed region left of the 3' SID,
// so matched bases pull the 3' trim boundary earlier (trim_pos decreases). Note the terminal-state
// groupings differ from the forward case -- here kInputEnd is scored as deletions (no read bases
// consumed) while kFoundPruned/kUnresolved are scored as substitutions that move the trim boundary.
ScoreTrim FoldBackwardExtension(const DFAClassifier::DFAResults& bw, s32 score, s32 trim_pos) {
  using enum DFAClassifier::DFAResult::State;
  if (!bw.last_match.has_value()) {
    const auto remaining_bases = bw.remaining_last_match_bases;
    return {score + (remaining_bases * scoring::kSubstitution), trim_pos - remaining_bases};
  }
  const auto& last = bw.last_match.value();
  trim_pos -= last.pos;
  const auto remaining_bases = bw.remaining_last_match_bases;
  switch (last.state) {
    case kFoundEnd:
      score += last.score;
      break;
    case kInputEnd:
      // Ran off the (likely truncated) end; the unmatched reference bases must be deletions.
      score += last.score + (remaining_bases * scoring::kDeletion);
      break;
    case kFoundPruned:
    case kUnresolved:
      score += last.score + (remaining_bases * scoring::kSubstitution);
      trim_pos -= remaining_bases;
      break;
    default:
      throw error::Error("Trim location ambiguous. Found state end {} ", magic_enum::enum_name(last.state));
  }
  return {score, trim_pos};
}

}  // namespace

DemuxAndTrimSimplex::DemuxAndTrimSimplex(const DemuxAndTrimParam& param,
                                         const vec<AdapterComponent>& components,
                                         DualSidMatchingPolicy dual_sid_matching_policy)
    : DemuxAndTrimSimplex{param, components, std::move(dual_sid_matching_policy), kIndependentStemMaxDrop} {
}

DemuxAndTrimSimplex::DemuxAndTrimSimplex(const DemuxAndTrimParam& param,
                                         const vec<AdapterComponent>& components,
                                         DualSidMatchingPolicy dual_sid_matching_policy,
                                         const s32 stem_max_drop)
    : _param{param},
      _sid_fwd_matcher{MatcherAt(components, kSidFwd)},
      _sid_rev_matcher{MatcherAt(components, kSidRev)},
      _dual_sid_matching_policy{std::move(dual_sid_matching_policy)},
      _nominal_seed_size_5p{static_cast<s32>(_sid_fwd_matcher->Lut()->Pool().front().sequence.size())},
      _nominal_seed_size_3p{static_cast<s32>(_sid_rev_matcher->Lut()->Pool().front().sequence.size())},
      // The OUTER anchors face the read ends: capping them at kExtendAnchorDfaMaxLen (with per-base
      // backoff in MakeExtendDfa) keeps the long DSL anchors (e.g. the 29-34bp I2X-S-DX anchors)
      // within the DFA powerset-state limit and only limits how far trimming reaches towards the read
      // end. The INNER stems face the insert, so they must be built at FULL length -- capping them
      // would leave the stem tail attached to the demuxed insert (the dual-index stem-trimming bug).
      // The per-region max_drop constants are defined above.
      _bw_extend_anchor_fwd_dfa{MakeExtendDfa(FixedAt(components, kAnchorFwd).Seq(),
                                              SearchDirection::kBackward,
                                              kBwExtendAnchorFwdMaxDrop,
                                              kExtendAnchorDfaMaxLen)},
      _fw_extend_stem_fwd_dfa{MakeExtendDfa(
          FixedAt(components, kStemFwd).Seq(), SearchDirection::kForward, stem_max_drop, kExtendStemDfaMaxLen)},
      _bw_extend_stem_rev_dfa{MakeExtendDfa(
          FixedAt(components, kStemRev).Seq(), SearchDirection::kBackward, stem_max_drop, kExtendStemDfaMaxLen)},
      _fw_extend_anchor_rev_dfa{MakeExtendDfa(FixedAt(components, kAnchorRev).Seq(),
                                              SearchDirection::kForward,
                                              kFwExtendAnchorRevMaxDrop,
                                              kExtendAnchorDfaMaxLen)},
      _fw_search_anchor_fwd_bitap(FixedAt(components, kAnchorFwd).Seq(), SearchDirection::kForward),
      _bw_search_anchor_rev_bitap(FixedAt(components, kAnchorRev).Seq(), SearchDirection::kBackward),
      _fw_partial_stem_fwd_bitap(FixedAt(components, kStemFwd).Seq(), SearchDirection::kForward),
      _bw_partial_stem_fwd_bitap(FixedAt(components, kStemFwd).Seq(), SearchDirection::kBackward),
      _fw_partial_stem_rev_bitap(FixedAt(components, kStemRev).Seq(), SearchDirection::kForward),
      _bw_partial_stem_rev_bitap(FixedAt(components, kStemRev).Seq(), SearchDirection::kBackward),
      _sid_5p_search_space_coefficients(sid_scoring::SidScoring::MakeSearchSpaceCoefficients(*_sid_fwd_matcher)),
      _sid_3p_search_space_coefficients(sid_scoring::SidScoring::MakeSearchSpaceCoefficients(*_sid_rev_matcher)) {
  if (_dual_sid_matching_policy.Enabled() &&
      (_dual_sid_matching_policy.sample_by_5p_sid.size() != _sid_fwd_matcher->Pool().size() ||
       _dual_sid_matching_policy.sample_by_3p_sid.size() != _sid_rev_matcher->Pool().size())) {
    throw error::Error("Dual-SID matching policy side-pool sizes do not match the simplex SID matchers");
  }

  ValidateMinReadLength(param, _sid_5p_search_space_coefficients, "5p SID");
  ValidateMinReadLength(param, _sid_3p_search_space_coefficients, "3p SID");

  Logging::Debug(
      "Running adapter design with min_score of {} ({} 5' / {} 3' sids). At {} bp the 5'/3' score "
      "thresholds are {}/{}; at {} bp they are {}/{}.",
      param.min_score,
      _sid_fwd_matcher->Pool().size(),
      _sid_rev_matcher->Pool().size(),
      param.min_read_len,
      param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(_sid_5p_search_space_coefficients,
                                                                        static_cast<u32>(param.min_read_len)),
      param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(_sid_3p_search_space_coefficients,
                                                                        static_cast<u32>(param.min_read_len)),
      kMaxReadLength,
      param.min_score +
          sid_scoring::SidScoring::CalcSearchSpacePenalty(_sid_5p_search_space_coefficients, kMaxReadLength),
      param.min_score +
          sid_scoring::SidScoring::CalcSearchSpacePenalty(_sid_3p_search_space_coefficients, kMaxReadLength));
}

/**
 * @brief Converts per-k unique k-mer counts into random match probabilities.
 *
 * @param unique_kmer_counts  Number of unique k-mers for each k in [min_k, max_k], indexed by (k - min_k).
 * @param min_k  Minimum k-mer size.
 * @param max_k  Maximum k-mer size.
 * @return Vector of probabilities indexed by (k − min_k).
 */
static vec<f64> KmerCountsToProb(const vec<f64>& unique_kmer_counts, const u8 min_k, const u8 max_k) {
  vec<f64> sid_kmer_prob(max_k - min_k + 1);
  for (u8 k = min_k; k <= max_k; ++k) {
    const auto unique_count = unique_kmer_counts[k - min_k];
    const auto kmer_count = std::ldexp(1.0, 2 * k);
    sid_kmer_prob[k - min_k] = unique_count / kmer_count;
    Logging::Debug(
        "For k-mer size {}, unique k-mers in SIDs: {}, total possible k-mers: {}, random match probability: {}",
        k,
        static_cast<u64>(unique_count),
        static_cast<u64>(kmer_count),
        sid_kmer_prob[k - min_k]);
  }
  return sid_kmer_prob;
}

vec<f64> sid_scoring::SidScoring::CalcSIDProbRandomMatch(const SeqMatcher& sid_matcher,
                                                         const u8 min_k,
                                                         const u8 max_k) {
  vec<f64> counts(max_k - min_k + 1);
  for (u8 k = min_k; k <= max_k; ++k) {
    counts[k - min_k] = static_cast<f64>(sid_matcher.Lut()->HashTables().at(k).size());
  }
  return KmerCountsToProb(counts, min_k, max_k);
}

vec<f64> sid_scoring::SidScoring::CalcSIDProbRandomMatch(const SeqMatcher& sid_matcher_a,
                                                         const SeqMatcher& sid_matcher_b,
                                                         const u8 min_k,
                                                         const u8 max_k) {
  vec<f64> counts(max_k - min_k + 1);
  for (u8 k = min_k; k <= max_k; ++k) {
    counts[k - min_k] = static_cast<f64>(sid_matcher_a.Lut()->HashTables().at(k).size() +
                                         sid_matcher_b.Lut()->HashTables().at(k).size());
  }
  return KmerCountsToProb(counts, min_k, max_k);
}

sid_scoring::SidScoring::SearchSpaceCoefficients sid_scoring::SidScoring::CalcSearchSpaceCoefficients(
    const vec<f64>& sid_kmer_prob, const u8 min_k, const u8 max_k) {
  f64 linear = 0.0;
  f64 constant = 0.0;
  for (u8 k = min_k; k <= max_k; ++k) {
    const auto p = sid_kmer_prob[k - min_k];
    linear += p;
    constant += p * (1.0 - k);
  }
  Logging::Debug("Precomputed search-space coefficients: linear: {}, constant: {}", linear, constant);
  return {linear, constant};
}

sid_scoring::SidScoring::SearchSpaceCoefficients sid_scoring::SidScoring::MakeSearchSpaceCoefficients(
    const SeqMatcher& sid_matcher) {
  const auto min_k = CalcSidMinK(sid_matcher);
  const auto max_k = CalcSidMaxK(sid_matcher);
  const auto sid_kmer_prob = CalcSIDProbRandomMatch(sid_matcher, min_k, max_k);
  Logging::Debug("SID pool k-range for scoring: [{}, {}]", min_k, max_k);
  return CalcSearchSpaceCoefficients(sid_kmer_prob, min_k, max_k);
}

sid_scoring::SidScoring::SearchSpaceCoefficients sid_scoring::SidScoring::AddSearchSpaceCoefficients(
    const SearchSpaceCoefficients& first, const SearchSpaceCoefficients& second) {
  return {.linear = first.linear + second.linear, .constant = first.constant + second.constant};
}

sid_scoring::SidScoring::SearchSpaceCoefficients sid_scoring::SidScoring::MakePairedSearchSpaceCoefficients(
    const SeqMatcher& sid_matcher_5p, const SeqMatcher& sid_matcher_3p) {
  const auto min_k = CalcSidMinK(sid_matcher_5p);
  const auto max_k = CalcSidMaxK(sid_matcher_5p);
  if (min_k != CalcSidMinK(sid_matcher_3p) || max_k != CalcSidMaxK(sid_matcher_3p)) {
    throw error::Error("Paired simplex SIDs must have the same length and maximum edit distance");
  }
  return CalcSearchSpaceCoefficients(
      CalcSIDProbRandomMatch(sid_matcher_5p, sid_matcher_3p, min_k, max_k), min_k, max_k);
}

s32 sid_scoring::SidScoring::CalcSearchSpacePenalty(const SearchSpaceCoefficients& coefficients, const u32 seq_len) {
  const auto expected_random_hits = (coefficients.linear * seq_len) + coefficients.constant;
  return static_cast<s32>(std::ceil(std::log2(expected_random_hits > 1.0 ? expected_random_hits : 1.0)));
}

void sid_scoring::SidScoring::ValidateSidBounds(const SeqMatcher& sid_matcher) {
  const auto sid_len = sid_matcher.GetFrontSeq().size();
  const auto sid_edist = sid_matcher.GetMaxDist();
  if (sid_len < sid_edist) {
    throw error::Error("SID length {} is shorter than edit distance {}", sid_len, sid_edist);
  }
}

u8 sid_scoring::SidScoring::CalcSidMinK(const SeqMatcher& sid_matcher) {
  ValidateSidBounds(sid_matcher);
  return static_cast<u8>(sid_matcher.GetFrontSeq().size() - sid_matcher.GetMaxDist());
}

u8 sid_scoring::SidScoring::CalcSidMaxK(const SeqMatcher& sid_matcher) {
  return static_cast<u8>(sid_matcher.GetFrontSeq().size() + sid_matcher.GetMaxDist());
}

DemuxAndTrimSimplex::SearchBounds DemuxAndTrimSimplex::MakeSearchBounds(const s32 length) const {
  return {
      .length = length,
      // Assignment must not depend on --min-trimmed-read-length. Keep only the non-empty insert
      // invariant here; the configured minimum is applied by the post-assignment read filters.
      .min_insert_len = 1,
      .min_score_5p = _param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(
                                             _sid_5p_search_space_coefficients, static_cast<u32>(length)),
      .min_score_3p = _param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(
                                             _sid_3p_search_space_coefficients, static_cast<u32>(length)),
  };
}

std::optional<DemuxAndTrimSimplex::ScoredMatch> DemuxAndTrimSimplex::Score5pMatch(const FixedReadRecord& record,
                                                                                  const MatchInfo& match,
                                                                                  const SearchBounds& bounds) const {
  if (match.SPos() > match.EPos() || match.EPos() > record.SeqLen()) {
    return std::nullopt;
  }
  auto dfa_results = ExtendCandidate(record, match, _bw_extend_anchor_fwd_dfa, _fw_extend_stem_fwd_dfa);
  const auto seed_score = match.LogLikelihoodScore(_nominal_seed_size_5p) + dfa_results.bw_dfa_results.best.score;
  const auto [score, trim_pos] =
      FoldForwardExtension(dfa_results.fw_dfa_results, seed_score, static_cast<s32>(match.EPos()));
  if (trim_pos < 0 || trim_pos > bounds.length - bounds.min_insert_len) {
    return std::nullopt;
  }
  return ScoredMatch{match, score, trim_pos};
}

std::optional<DemuxAndTrimSimplex::ScoredMatch> DemuxAndTrimSimplex::Score3pMatch(const FixedReadRecord& record,
                                                                                  const MatchInfo& match,
                                                                                  const SearchBounds& bounds) const {
  if (match.SPos() > match.EPos() || match.EPos() > record.SeqLen()) {
    return std::nullopt;
  }
  auto dfa_results = ExtendCandidate(record, match, _bw_extend_stem_rev_dfa, _fw_extend_anchor_rev_dfa);
  const auto seed_score = match.LogLikelihoodScore(_nominal_seed_size_3p) + dfa_results.fw_dfa_results.best.score;
  const auto [score, trim_pos] =
      FoldBackwardExtension(dfa_results.bw_dfa_results, seed_score, static_cast<s32>(match.SPos()));
  if (trim_pos < bounds.min_insert_len || trim_pos > bounds.length) {
    return std::nullopt;
  }
  return ScoredMatch{match, score, trim_pos};
}

std::optional<DemuxAndTrimSimplex::ScoredMatch> DemuxAndTrimSimplex::FindBest5p(const FixedReadRecord& record,
                                                                                const SearchBounds& bounds,
                                                                                const s32 search_limit) const {
  std::optional<ScoredMatch> best;
  auto search_pos = 0;
  vec<MatchInfo> matches;
  while (search_pos < search_limit) {
    _sid_fwd_matcher->FindNextBarcode(
        ReadEnd::k5p, search_pos, record.TwoBitsSeq(), static_cast<size_t>(search_limit), matches);
    if (matches.empty() || matches.front().SPos() >= static_cast<u32>(search_limit)) {
      break;
    }
    search_pos = static_cast<s32>(matches.front().SPos()) + 1;
    for (const auto& match : matches) {
      auto candidate = Score5pMatch(record, match, bounds);
      if (candidate && (!best || candidate->score > best->score)) {
        best = candidate;
      }
    }
  }
  return best;
}

std::optional<DemuxAndTrimSimplex::ScoredMatch> DemuxAndTrimSimplex::FindBest3p(const FixedReadRecord& record,
                                                                                const SearchBounds& bounds,
                                                                                const s32 search_limit) const {
  std::optional<ScoredMatch> best;
  auto search_pos = bounds.length;
  vec<MatchInfo> matches;
  while (search_pos > search_limit) {
    _sid_rev_matcher->FindNextBarcode(
        ReadEnd::k3p, search_pos, record.TwoBitsSeq(), static_cast<size_t>(bounds.length), matches);
    if (matches.empty() || matches.front().EPos() == 0 || matches.front().EPos() > static_cast<u32>(bounds.length)) {
      break;
    }
    search_pos = static_cast<s32>(matches.front().EPos()) - 1;
    for (const auto& match : matches) {
      auto candidate = Score3pMatch(record, match, bounds);
      if (candidate && (!best || candidate->score > best->score)) {
        best = candidate;
      }
    }
  }
  return best;
}

std::optional<LociRange> DemuxAndTrimSimplex::FindTerminalStem(const FixedReadRecord& record, const bool is_5p) const {
  const auto& forward = is_5p ? _fw_partial_stem_fwd_bitap : _fw_partial_stem_rev_bitap;
  const auto& backward = is_5p ? _bw_partial_stem_fwd_bitap : _bw_partial_stem_rev_bitap;
  const auto& anchor = is_5p ? _fw_search_anchor_fwd_bitap : _bw_search_anchor_rev_bitap;
  const auto unresolved_sid_length =
      is_5p ? _sid_fwd_matcher->GetFrontSeq().size() : _sid_rev_matcher->GetFrontSeq().size();
  const auto read_end = static_cast<s32>(record.SeqLen()) - 1;
  const auto search_begin = is_5p ? 0 : std::max(0, read_end - kPartialStemEndpointWindow + 1);
  const auto search_end = is_5p ? std::min(read_end, kPartialStemEndpointWindow - 1) : read_end;
  const auto scan_step =
      Bitap<kPartialStemMaxEditDistance>::kQueryWindowSize - forward.GetQueryLength() - kPartialStemMaxEditDistance;

  std::optional<PartialStemCandidate> best;

  // Bitap operates on one 64bp word. Overlap each window by the full query plus edit band so a
  // valid 33/34bp alignment is never split across windows. The directional searches only propose
  // endpoints: validate the complete resulting interval below before it can affect trimming.
  for (s32 window_begin = search_begin; window_begin <= search_end; window_begin += scan_step) {
    const auto window_end =
        std::min(window_begin + Bitap<kPartialStemMaxEditDistance>::kQueryWindowSize - 1, search_end);
    const auto candidate = FindTerminalStemCandidate(record, is_5p, forward, backward, window_begin, window_end);
    if (!candidate.has_value()) {
      continue;
    }

    const auto outer_context =
        is_5p ? std::string_view(record.Seq(), candidate->range.spos)
              : std::string_view(record.Seq() + candidate->range.epos, record.SeqLen() - candidate->range.epos);
    if (!HasPlausibleOuterContext(anchor.Query(), outer_context, unresolved_sid_length, is_5p)) {
      continue;
    }

    // Do not sum the two directional Bitap distances: they can describe the same edits. The exact
    // interval distance is both the acceptance criterion and the only comparable candidate score.
    if (!best.has_value() || IsBetterPartialStemCandidate(*candidate, *best, is_5p)) {
      best = candidate;
    }
  }
  return best.has_value() ? std::optional<LociRange>{best->range} : std::nullopt;
}

std::optional<LociRange> DemuxAndTrimSimplex::FindTerminalPartialStem(const FixedReadRecord& record,
                                                                      const bool is_5p) const {
  const auto stem = is_5p ? _fw_partial_stem_fwd_bitap.Query() : _fw_partial_stem_rev_bitap.Query();

  // A cutoff 5' stem must start at the physical 5' read boundary, so only its insert-facing
  // suffix can remain. A cutoff 3' stem must end at the physical 3' read boundary, so only its
  // insert-facing prefix can remain. Do not scan the endpoint window here: that could trim a
  // coincidental stem-like substring inside the insert.
  // A likelihood-aware stem score calibrated against negative controls can eventually replace this
  // edit-distance acceptance rule and give terminal partial-stem trimming an explicit false-positive rate.
  const auto read = std::string_view(record.Seq(), record.SeqLen());
  const auto candidate = FindTerminalPartialStemCandidate(stem, read, is_5p);
  return candidate.has_value() ? std::optional<LociRange>{candidate->range} : std::nullopt;
}

void DemuxAndTrimSimplex::CompletePartialStem(const FixedReadRecord& record,
                                              const SearchBounds& bounds,
                                              TrimInfoSimplex& result) const {
  if (!_dual_sid_matching_policy.Enabled() || !result.sid.has_value() ||
      result.sid_5p.has_value() == result.sid_3p.has_value()) {
    return;
  }

  if (result.sid_5p.has_value()) {
    const auto stem_3p =
        FindTerminalStem(record, false).or_else([this, &record] { return FindTerminalPartialStem(record, false); });
    if (stem_3p && result.insert.spos + static_cast<u32>(bounds.min_insert_len) <= stem_3p->spos) {
      result.insert.epos = stem_3p->spos;
    }
  } else {
    const auto stem_5p =
        FindTerminalStem(record, true).or_else([this, &record] { return FindTerminalPartialStem(record, true); });
    if (stem_5p && stem_5p->epos + static_cast<u32>(bounds.min_insert_len) <= result.insert.epos) {
      result.insert.spos = stem_5p->epos;
    }
  }
}

void DemuxAndTrimSimplex::ReconstructIndependentResult(const SearchBounds& bounds,
                                                       const std::optional<ScoredMatch>& best_5p,
                                                       const std::optional<ScoredMatch>& best_3p,
                                                       TrimInfoSimplex& result) const {
  const bool has_5p = best_5p && best_5p->score >= bounds.min_score_5p;
  const bool has_3p = best_3p && best_3p->score >= bounds.min_score_3p;
  if (has_5p) {
    result.sid_5p = best_5p->match.BarcodeId();
    result.sid_5p_edist = best_5p->match.EDist();
    result.score_5p = best_5p->score;
  }
  if (has_3p) {
    result.sid_3p = best_3p->match.BarcodeId();
    result.sid_3p_edist = best_3p->match.EDist();
    result.score_3p = best_3p->score;
  }

  if (has_5p && has_3p) {
    if (best_5p->trim_pos < best_3p->trim_pos) {
      const auto sample_id = PairedSampleId(best_5p->match.BarcodeId(), best_3p->match.BarcodeId());
      if (sample_id.has_value()) {
        result.sid = sample_id;
        result.sid_5p = sample_id;
        result.sid_3p = sample_id;
        result.insert.spos = static_cast<u32>(best_5p->trim_pos);
        result.insert.epos = static_cast<u32>(best_3p->trim_pos);
      }
    }
    return;
  }

  if (has_3p) {
    const auto sample_id = SampleIdFor3pPartialMatch(best_3p->match.BarcodeId());
    if (sample_id.has_value()) {
      result.sid = sample_id;
      result.sid_3p = sample_id;
      result.insert.epos = static_cast<u32>(best_3p->trim_pos);
    }
  } else if (has_5p) {
    const auto sample_id = SampleIdFor5pPartialMatch(best_5p->match.BarcodeId());
    if (sample_id.has_value()) {
      result.sid = sample_id;
      result.sid_5p = sample_id;
      result.insert.spos = static_cast<u32>(best_5p->trim_pos);
      result.insert.epos = static_cast<u32>(bounds.length);
    }
  }
}

const TrimInfoSimplex& DemuxAndTrimSimplex::operator()(FixedReadRecord& record) const {
  auto& result = record.trim_info_simplex;
  result.Clear();
  const auto bounds = MakeSearchBounds(static_cast<s32>(record.SeqLen()));
  result.insert.epos = static_cast<u32>(bounds.length);
  const auto best_5p = FindBest5p(record, bounds, bounds.length - bounds.min_insert_len);
  const auto best_3p = FindBest3p(record, bounds, bounds.min_insert_len);
  ReconstructIndependentResult(bounds, best_5p, best_3p, result);
  CompletePartialStem(record, bounds, result);
  return result;
}

DemuxAndTrimSimplex::PairedDFAResults DemuxAndTrimSimplex::ExtendCandidate(const FixedReadRecord& record,
                                                                           const MatchInfo& match_info,
                                                                           const DFAClassifier& bw_dfa,
                                                                           const DFAClassifier& fw_dfa) {
  const auto search_domain_bw = std::string_view(record.Seq(), match_info.SPos());
  const auto search_domain_fw = std::string_view(record.Seq() + match_info.EPos(), record.SeqLen() - match_info.EPos());
  return {bw_dfa.Classify(search_domain_bw), fw_dfa.Classify(search_domain_fw)};
}

void DemuxAndTrimSimplex::ValidateMinReadLength(const DemuxAndTrimParam& param,
                                                const sid_scoring::SidScoring::SearchSpaceCoefficients& coefficients,
                                                const std::string_view scoring_mode) {
  const auto threshold = param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(
                                               coefficients, static_cast<u32>(param.min_read_len));
  const auto min_match_length = static_cast<u32>((threshold + scoring::kMatch - 1) / scoring::kMatch);
  if (param.min_read_len < min_match_length) {
    throw error::Error("The {} requires at least {} searchable bases at --min-score {}, but --min-read-length is {}",
                       scoring_mode,
                       min_match_length,
                       param.min_score,
                       param.min_read_len);
  }
}

std::optional<u32> DemuxAndTrimSimplex::PairedSampleId(const u32 sid_5p_id, const u32 sid_3p_id) const {
  if (!_dual_sid_matching_policy.Enabled()) {
    return sid_5p_id == sid_3p_id ? std::optional<u32>{sid_5p_id} : std::nullopt;
  }
  const auto it = _dual_sid_matching_policy.sample_by_sid_pair.find(DualSidPairKey(sid_5p_id, sid_3p_id));
  if (it == _dual_sid_matching_policy.sample_by_sid_pair.end()) {
    return std::nullopt;
  }
  return it->second;
}

std::optional<u32> DemuxAndTrimSimplex::SampleIdFor5pPartialMatch(const u32 sid_5p_id) const {
  if (!_dual_sid_matching_policy.Enabled()) {
    return sid_5p_id;
  }
  if (sid_5p_id >= _dual_sid_matching_policy.sample_by_5p_sid.size()) {
    return std::nullopt;
  }
  const auto sample_id = _dual_sid_matching_policy.sample_by_5p_sid.at(sid_5p_id);
  if (!sample_id.has_value() || *sample_id >= _dual_sid_matching_policy.partial_matches.size()) {
    return std::nullopt;
  }
  return _dual_sid_matching_policy.partial_matches.at(*sample_id).allow_5p_partial_match ? sample_id : std::nullopt;
}

std::optional<u32> DemuxAndTrimSimplex::SampleIdFor3pPartialMatch(const u32 sid_3p_id) const {
  if (!_dual_sid_matching_policy.Enabled()) {
    return sid_3p_id;
  }
  if (sid_3p_id >= _dual_sid_matching_policy.sample_by_3p_sid.size()) {
    return std::nullopt;
  }
  const auto sample_id = _dual_sid_matching_policy.sample_by_3p_sid.at(sid_3p_id);
  if (!sample_id.has_value() || *sample_id >= _dual_sid_matching_policy.partial_matches.size()) {
    return std::nullopt;
  }
  return _dual_sid_matching_policy.partial_matches.at(*sample_id).allow_3p_partial_match ? sample_id : std::nullopt;
}

}  // namespace xoos::demux
