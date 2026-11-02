#pragma once

#include <limits>
#include <optional>

#include "adapters/simplex-10x/lut-bundle-simplex-10x.h"
#include "adapters/simplex-10x/sid-fast-match.h"
#include "adapters/simplex/demux-and-trim-simplex.h"
#include "adapters/simplex/trim-info-simplex.h"
#include "core/demux-and-trim-pipeline.h"
#include "io/read-record.h"
#include "sequence/matcher/bitap.h"
#include "sequence/matcher/dfa-classifier.h"

namespace xoos::demux {

/**
 * @class DemuxAndTrimSimplex10x
 * @brief Simplex-10x trimmer with hash-table fast-path and LUT+DFA+Bitap slow-path.
 *
 * Read structure: [fixed_1 (runway)][sid_2 (SID 12 bp)][fixed_3 (stem/p_read_1t 22 bp)][insert...]
 *
 * Only the runway and SID are trimmed. The trim position is set after the SID —
 * the stem is left in the insert.
 *
 * Fast-path (~86% of reads): exact runway match at offset 0 or 1, followed by
 * a 2-bit encoded hash-table lookup of (SID + 4bp stem prefix). Handles exact
 * SID matches and 1-substitution variants. The 16bp key provides ~1/4^16
 * specificity per position.
 *
 * Slow-path (fallback): LUT seed scan + backward DFA extension into fixed_1
 * (runway) + Bitap<4> stem confirmation + Karlin-Altschul log-likelihood
 * scoring with search-space penalty. Per-SID best-score aggregation prevents
 * variable-length k-mer matches from distorting scores.
 */
class DemuxAndTrimSimplex10x {
 public:
  DemuxAndTrimSimplex10x(const DemuxAndTrimParam& param, const LutBundleSimplex10x& lut_bundle);

  /**
   * @brief Initialize the fast-path using the full runway (e.g. "CAACAA").
   *
   * TryMatch checks for the truncated runway (runway[1:]) at offset 0
   * first, then the full runway at offset 0.
   *
   * @param runway    Full runway sequence (e.g. "CAACAA").
   * @param sid_pool  SID barcode pool.
   * @param stem      Full stem sequence.
   */
  void InitFastPath(std::string_view runway, const BarcodePool& sid_pool, std::string_view stem);

  /**
   * @brief Identifies the SID barcode and computes 5' trim boundaries.
   * @param record The read record to demultiplex.
   * @return TrimInfoSimplex with SID assignment and insert coordinates.
   */
  TrimInfoSimplex operator()(const FixedReadRecord& record) const;

 private:
  /** @brief Per-SID best-score tracking used during the LUT scan phase. */
  struct SidCandidate {
    s32 score = std::numeric_limits<s32>::min();
    MatchInfo match;
    s32 trim_pos = 0;
    s32 hit_count = 0;
  };

  /**
   * @brief Try the fast-path if it is initialized.
   * @return The match result, or nullopt when the fast-path is not available
   *         or the read does not match.
   */
  std::optional<SidFastMatch::Result> TryFastPath(const char* seq, size_t length) const;

  /** @brief Build a 5'-only TrimInfoSimplex from a SID match result and read length. */
  static TrimInfoSimplex MakeTrimResult(const SidFastMatch::Result& match, u32 read_len);
  static TrimInfoSimplex MakeTrimResult(const SidCandidate& candidate, u32 read_len);

  /** @brief Score a single SID seed match (steps 2a–2d) and update per_sid tracking. */
  void ScoreSidCandidate(const MatchInfo& sid_match,
                         const FixedReadRecord& record,
                         std::vector<SidCandidate>& per_sid) const;

  /**
   * @brief Compute log-likelihood seed score with the match length capped to
   * prevent flanking-base absorption from inflating the score.
   */
  static s32 CappedLogLikelihoodScore(const MatchInfo& match, s32 capped_len, s32 gt_len);

  const DemuxAndTrimParam& _param;
  const LutBundleSimplex10x& _lut_bundle;

  const s32 _nominal_seed_size;

  // Backward DFA for runway (fixed_1) extension scoring.
  const DFAClassifier _bw_extend_1_dfa;

  // Backward Bitap<4> for fixed_3 (stem/p_read_1t). A single backward pass
  // returns the match start position directly, handling indels correctly.
  const Bitap<4> _fixed_3_bitap;
  const s32 _fixed_3_len;

  // Maximum 5' trim position for the LUT scan, derived from the adapter
  // structure: fixed_1 + SID + margin for linker bases and sequencing errors.
  const s32 _max_adapter_span;

  const u8 _sid_min_k;
  const u8 _sid_max_k;
  const std::vector<f64> _sid_kmer_prob;
  const sid_scoring::SidScoring::SearchSpaceCoefficients _search_space_coefficients;

  // Fast-path: single matcher using the full runway (e.g. CAACAA).
  // TryMatch checks runway[1:] at offset 0 first, then the full runway.
  std::optional<SidFastMatch> _fast_match;
};
}  // namespace xoos::demux
