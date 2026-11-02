#include "demux-and-trim-simplex-10x.h"

#include "xoos/log/logging.h"

namespace xoos::demux {

// Margin added to the nominal adapter span (fixed_1 + SID + edit distance)
// when computing the LUT scan cap. Breakdown:
//   - Linker bases between SID and fixed_3:  0-2 bp
//   - Possible leading base before fixed_1:  0-1 bp
//   - Sequencing-error tolerance (indels,    ~9 bp
//     homopolymer shifts, adapter damage)
// Total: 12 bp. This is deliberately generous — a too-tight cap silently
// drops valid hits, while a too-loose cap only costs LUT scan cycles.
constexpr s32 kScanMargin = 12;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
// Initializes the scoring infrastructure and matchers for the simplex-10x
// adapter. The adapter structure is 5'-only:
//
//   [fixed_1 (runway)] [sid_2 (SID barcode, 12 bp)] [fixed_3 (stem, 22 bp)] [insert...]
//
// The runway, SID, and any linker bases are trimmed. The trim position is
// set after the SID — the stem is left in the insert.
//
// The slow-path reuses the LUT seed scan + backward DFA extension +
// Bitap<4> stem confirmation + Karlin-Altschul scoring approach from the
// simplex adapter (see demux-and-trim-simplex.cpp). The fast-path
// (SidFastMatch) is initialized separately via InitFastPath().
DemuxAndTrimSimplex10x::DemuxAndTrimSimplex10x(const DemuxAndTrimParam& param, const LutBundleSimplex10x& lut_bundle)
    : _param{param},
      _lut_bundle{lut_bundle},
      _nominal_seed_size{static_cast<s32>(lut_bundle.sid_2_matcher.Lut()->Pool().front().sequence.size())},
      _bw_extend_1_dfa{_lut_bundle.fixed_1.Seq(), SearchDirection::kBackward},
      _fixed_3_bitap{_lut_bundle.fixed_3.Seq(), SearchDirection::kBackward},
      _fixed_3_len{static_cast<s32>(_lut_bundle.fixed_3.Size())},
      // LUT scan cap: fixed_1 + SID + slack for linker (0-2 bp), leading base,
      // and edit-distance-induced position shifts in the SID match.
      _max_adapter_span{static_cast<s32>(_lut_bundle.fixed_1.Size()) + _nominal_seed_size +
                        static_cast<s32>(_lut_bundle.sid_2_matcher.GetMaxDist()) + kScanMargin},
      _sid_min_k(sid_scoring::SidScoring::CalcSidMinK(lut_bundle.sid_2_matcher)),
      _sid_max_k(
          static_cast<u8>(_lut_bundle.sid_2_matcher.GetFrontSeq().size() + _lut_bundle.sid_2_matcher.GetMaxDist())),
      _sid_kmer_prob(sid_scoring::SidScoring::CalcSIDProbRandomMatch(lut_bundle.sid_2_matcher, _sid_min_k, _sid_max_k)),
      _search_space_coefficients(
          sid_scoring::SidScoring::CalcSearchSpaceCoefficients(_sid_kmer_prob, _sid_min_k, _sid_max_k)) {
  // Validate that the minimum read length is long enough to produce a
  // meaningful score above the threshold (same check as simplex adapter).
  const auto min_length_score_threshold =
      param.min_score +
      sid_scoring::SidScoring::CalcSearchSpacePenalty(_search_space_coefficients, static_cast<u32>(param.min_read_len));

  const auto min_match_length = static_cast<u32>((min_length_score_threshold + scoring::kMatch - 1) / scoring::kMatch);

  if (param.min_read_len < min_match_length) {
    throw error::Error(
        "The minimum adapter length is {}, but --min_read_length is set to {}. "
        "Set --min_read_length to a value >= {} for this adapter.",
        min_match_length,
        param.min_read_len,
        min_match_length);
  }

  Logging::Debug(
      "Running simplex-10x adapter (LUT+DFA SID, Bitap<4> stem scoring) with min_score of {} and {} sids. "
      "Score threshold of {} for sequences {} bp long and {} for sequences {} bp long.",
      param.min_score,
      _lut_bundle.sid_2_matcher.Pool().size(),
      min_length_score_threshold,
      param.min_read_len,
      param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(_search_space_coefficients, kMaxReadLength),
      kMaxReadLength);
}

// ---------------------------------------------------------------------------
// Search-space penalty (Karlin-Altschul framework)
// ---------------------------------------------------------------------------
// The search-space penalty (computed via SidScoring::CalcSearchSpacePenalty)
// accounts for the probability of a random k-mer match in a read of a given
// length. This prevents long reads (which have more chances for spurious LUT
// hits) from producing inflated scores. The approach is identical to the
// simplex adapter (see demux-and-trim-simplex.cpp).

/**
 * Computes a log-likelihood seed score with the match length capped at the
 * nominal SID size. The LUT's variable-length k-mer matching can return
 * windows longer than the SID when flanking bases (runway or stem) are
 * absorbed into the match. Capping prevents this absorption from inflating
 * the score. The decomposition into matches/substitutions/insertions/deletions
 * follows the same logic as MatchInfo::LogLikelihoodScore.
 *
 * Note: EDist() is used uncapped because the LUT only absorbs flanking bases
 * that are exact matches — the edit distance reflects errors within the SID
 * itself, not in the absorbed flanking region.
 */
s32 DemuxAndTrimSimplex10x::CappedLogLikelihoodScore(const MatchInfo& match, const s32 capped_len, const s32 gt_len) {
  const s32 net_indels = capped_len - gt_len;
  const s32 abs_net = net_indels < 0 ? -net_indels : net_indels;

  const s32 num_insertions = net_indels > 0 ? net_indels : 0;
  const s32 num_deletions = net_indels < 0 ? -net_indels : 0;
  const s32 num_substitutions = static_cast<s32>(match.EDist()) - abs_net;
  const s32 num_matches = gt_len - num_substitutions - num_deletions;

  return (num_matches * scoring::kMatch) + (num_substitutions * scoring::kSubstitution) +
         (num_insertions * scoring::kInsertion) + (num_deletions * scoring::kDeletion);
}

// ---------------------------------------------------------------------------
// ScoreSidCandidate — inner loop body for the LUT scan
// ---------------------------------------------------------------------------
// Scores a single SID seed match by combining:
//   (a) Backward DFA extension into fixed_1 (runway).
//   (b) Capped seed log-likelihood score.
//   (c) Bitap<4> stem confirmation and trim position determination.
// The trim position is set to the start of fixed_3 (stem / p_read_1t) when
// found, ensuring the full SID and any linker bases are removed from the
// insert. Falls back to the SID-based estimate when the stem is not found.
// Updates the per-SID best-score tracker if this hit improves the score.
void DemuxAndTrimSimplex10x::ScoreSidCandidate(const MatchInfo& sid_match,
                                               const FixedReadRecord& record,
                                               std::vector<SidCandidate>& per_sid) const {
  // --- Step 2a: Backward DFA extension into fixed_1 (runway) ---
  //
  // The DFA scores how well the sequence upstream of the SID matches the
  // expected runway pattern (e.g., "CAACAA"). The domain is Read[0:SPos].
  //
  // Normalization: the LUT may return match windows longer than the
  // nominal SID size when flanking bases are absorbed into the k-mer.
  // Use max(SPos, EPos - nominal_sid_size) as the domain boundary so
  // that over-long matches don't shrink the runway domain and distort
  // the DFA score.
  const auto normalized_spos =
      std::max(static_cast<s32>(sid_match.SPos()), static_cast<s32>(sid_match.EPos()) - _nominal_seed_size);
  const auto search_domain_bw = std::string_view(record.Seq(), static_cast<size_t>(normalized_spos));
  const auto bw_dfa_results = _bw_extend_1_dfa.Classify(search_domain_bw);

  // --- Step 2b: Seed score with length capping ---
  //
  // Cap the seed match length at the nominal SID size for scoring. When
  // the LUT returns a window longer than the SID, the extra bases are
  // flanking sequence and should not inflate the seed score. Uses the
  // same log-likelihood decomposition as MatchInfo::LogLikelihoodScore
  // (matches * kMatch + subs * kSubstitution + ...).
  const auto capped_len = std::min(static_cast<s32>(sid_match.Length()), _nominal_seed_size);
  const auto seed_score = CappedLogLikelihoodScore(sid_match, capped_len, _nominal_seed_size);

  auto score = seed_score + bw_dfa_results.best.score;

  // --- Step 2c: Bitap<4> stem confirmation and trim position ---
  //
  // Search for fixed_3 (p_read_1t) starting right after the SID match end,
  // tolerating up to 4 edit distance. A single backward Bitap pass returns
  // the match start position directly, handling indels correctly without
  // needing a second pass. The trim point is the start of fixed_3 —
  // everything before it (runway, SID, and any linker bases) is trimmed.
  //
  // kStemSearchMargin: extra bases beyond fixed_3 length for the search
  // window right bound. Covers:
  //   - Linker bases between SID and fixed_3:  0-2 bp
  //   - SID-fixed_3 overlap from LUT k>12:     0-2 bp
  //   - Indel slack (Bitap ed=4):              0-4 bp
  // Total: up to 8 bp; rounded to 10 for safety.
  constexpr s32 kStemSearchMargin = 10;
  const auto search_start = static_cast<s32>(sid_match.EPos());
  const auto search_end =
      std::min(search_start + _fixed_3_len + kStemSearchMargin, static_cast<s32>(record.SeqLen() - 1));

  s32 stem_start = -1;
  if (search_start <= search_end) {
    const auto seq = std::string_view(record.Seq(), record.SeqLen());
    const s32 match_start = _fixed_3_bitap.Find(seq, search_start, search_end);
    if (match_start != -1) {
      score += _fixed_3_len * scoring::kMatch;
      stem_start = match_start;
    } else {
      score += _fixed_3_len * scoring::kSubstitution;
    }
  } else {
    score += _fixed_3_len * scoring::kSubstitution;
  }

  // Trim position: use the stem (fixed_3) start when found, otherwise
  // fall back to the SID-based estimate. The max() handles two LUT edge
  // cases:
  //   - EPos < nominal end: the LUT matched a shorter k-mer (k < 12),
  //     so EPos falls short of the true SID boundary.
  //   - EPos > nominal end: the LUT matched a longer k-mer (k > 12),
  //     absorbing linker or fixed_3 bases. normalized_spos + nominal_seed_size
  //     gives the correct boundary.
  // When Bitap finds the stem, stem_start supersedes both — it marks the
  // true fixed_3 start regardless of LUT window size. The sid_based_trim
  // fallback is only used when Bitap fails to find fixed_3 (e.g., heavily
  // damaged adapter), in which case trimming at the SID boundary is the
  // safest conservative choice.
  const auto sid_based_trim = std::max(static_cast<s32>(sid_match.EPos()), normalized_spos + _nominal_seed_size);
  const auto trim_pos_5p = (stem_start >= sid_based_trim) ? stem_start : sid_based_trim;

  // --- Step 2d: Update per-SID best score and hit count ---
  auto& candidate = per_sid[sid_match.BarcodeId()];
  ++candidate.hit_count;
  if (score > candidate.score) {
    candidate.score = score;
    candidate.match = sid_match;
    candidate.trim_pos = trim_pos_5p;
  }
}

// ---------------------------------------------------------------------------
// Main demux-and-trim operator
// ---------------------------------------------------------------------------
// Processes a single read to identify the SID barcode and determine the 5'
// trim position. The algorithm has four steps:
//
//   1. Compute search boundaries and the read-length-adjusted score threshold.
//   2. Scan the 5' region with the LUT to find SID seed matches. For each hit,
//      ScoreSidCandidate() computes the composite score and updates per-SID
//      tracking.
//   3. Select the globally best SID from the per-SID bests, using hit count
//      as a tiebreaker.
//   4. Assign the result if the best score passes the threshold.
//
// Key differences from the simplex adapter (demux-and-trim-simplex.cpp):
//   - 5'-only: no 3' SID search, no paired-SID logic.
//   - Stem (fixed_3) determines the trim position — the trim point is the
//     start of fixed_3, leaving the stem in the insert.
//   - Per-SID aggregation with DFA domain normalization and seed length
//     capping to handle the LUT's variable-length k-mer matching artifacts.
//
void DemuxAndTrimSimplex10x::InitFastPath(const std::string_view runway,
                                          const BarcodePool& sid_pool,
                                          const std::string_view stem) {
  _fast_match = SidFastMatch::Build(runway, sid_pool, stem);
}

TrimInfoSimplex DemuxAndTrimSimplex10x::MakeTrimResult(const SidFastMatch::Result& match, const u32 read_len) {
  TrimInfoSimplex result;
  result.sid = match.sid_id;
  result.sid_5p = match.sid_id;
  result.sid_5p_edist = match.edist;
  result.sid_3p = std::nullopt;
  result.sid_3p_edist = std::nullopt;
  result.insert = LociRange{match.insert_start, read_len};
  return result;
}

TrimInfoSimplex DemuxAndTrimSimplex10x::MakeTrimResult(const SidCandidate& candidate, const u32 read_len) {
  return MakeTrimResult({candidate.match.BarcodeId(), candidate.match.EDist(), static_cast<u32>(candidate.trim_pos)},
                        read_len);
}

std::optional<SidFastMatch::Result> DemuxAndTrimSimplex10x::TryFastPath(const char* const seq,
                                                                        const size_t length) const {
  if (_fast_match.has_value()) {
    return _fast_match->TryMatch(seq, length);
  }
  return std::nullopt;
}

TrimInfoSimplex DemuxAndTrimSimplex10x::operator()(const FixedReadRecord& record) const {
  const auto length = record.SeqLen();

  // Sanity check: the caller (Demux::operator()) already filters by
  // min_read_len, but guard here defensively.
  if (length < _param.min_read_len) {
    return {};
  }

  // Fast-path: hash-table SID lookup after exact runway match at offset 0.
  // The 16bp key (12bp SID + 4bp stem prefix) provides ~1/4^16 specificity
  // per position, making false positives negligible. No Bitap stem
  // confirmation needed — the stem prefix in the key is sufficient.
  auto fast_match_result = TryFastPath(record.Seq(), length);
  if (fast_match_result.has_value()) {
    return MakeTrimResult(*fast_match_result, length);
  }

  // --- Slow-path: LUT seed scan + backward DFA + Bitap stem confirmation ---

  // --- Step 1: Compute search boundaries and score threshold ---
  //
  // The LUT scan is bounded by _max_adapter_span (fixed_1 + SID + margin)
  // to avoid scanning deep into the insert.
  //
  // min_score: base threshold + search-space penalty for this read length.
  const auto min_score =
      _param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(_search_space_coefficients, length);

  TrimInfoSimplex result;
  result.insert.epos = length;

  // Per-SID best-score tracking. The LUT returns variable-length k-mer matches
  // that may absorb flanking sequence (runway or stem) into the SID match window,
  // shifting the match position and distorting downstream DFA scores. Tracking the
  // best score per SID ensures that a spurious cross-SID hit at a lucky position
  // can't beat the true SID's best evidence. Hit count serves as a tiebreaker:
  // the correct SID typically produces more seed hits than a coincidental cross-match.
  const auto pool_size = _lut_bundle.sid_2_matcher.Pool().size();
  std::vector<SidCandidate> per_sid(pool_size);

  // --- Step 2: Scan the 5' region for SID barcode seed matches ---
  //
  // The LUT hash table (SeqMatcher::FindNextBarcode) scans the read from left
  // to right, returning all SID candidates at each position. For each seed hit
  // ScoreSidCandidate() computes a composite score from three components:
  //   (a) Seed score: log-likelihood of the SID k-mer match.
  //   (b) Runway score: backward DFA extension into fixed_1.
  //   (c) Stem score: Bitap<4> confirmation of fixed_3 (scoring only, no trim).
  {
    auto search_pos = 0;
    std::vector<MatchInfo> sid_matches;

    const auto lut_scan_end = std::min(_max_adapter_span, static_cast<s32>(length));
    while (search_pos < lut_scan_end) {
      // Greedy search for the next barcode seed match. Returns all SIDs that
      // match at the same return position (EPos for 5' reads).
      _lut_bundle.sid_2_matcher.FindNextBarcode(
          ReadEnd::k5p, search_pos, record.TwoBitsSeq(), static_cast<size_t>(lut_scan_end), sid_matches);
      if (sid_matches.empty()) {
        break;
      }
      search_pos = static_cast<s32>(sid_matches[0].SPos()) + 1;
      for (const auto& sid_match : sid_matches) {
        ScoreSidCandidate(sid_match, record, per_sid);
      }
    }
  }

  // --- Step 3: Select the globally best SID from per-SID bests ---
  //
  // Primary criterion: highest composite score. Tiebreaker: more LUT hits
  // indicate stronger evidence (the correct SID typically produces multiple
  // seed matches at nearby positions while a spurious cross-match produces
  // few).
  auto best_score = std::numeric_limits<s32>::min();
  s32 best_hits = 0;
  const SidCandidate* best_candidate = nullptr;

  for (const auto& candidate : per_sid) {
    if (candidate.score > best_score || (candidate.score == best_score && candidate.hit_count > best_hits)) {
      best_score = candidate.score;
      best_hits = candidate.hit_count;
      best_candidate = &candidate;
    }
  }

  // --- Step 4: Assign the best candidate if it passes the score threshold ---
  //
  // Simplex-10x is 5'-only: sid_3p is always nullopt and insert.epos is the
  // full read length.
  if (best_candidate != nullptr && best_score >= min_score && best_candidate->match.HasMatches()) {
    return MakeTrimResult(*best_candidate, length);
  }

  return result;
}

}  // namespace xoos::demux
