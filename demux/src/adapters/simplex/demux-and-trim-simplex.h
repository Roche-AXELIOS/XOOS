#pragma once

#include <optional>
#include <string_view>

#include <xoos/types/vec.h>

#include "adapter-design/adapter-design-bundle.h"
#include "core/demux-and-trim-pipeline.h"
#include "core/dual-sid-mode.h"
#include "io/read-record.h"
#include "sequence/matcher/bitap.h"
#include "sequence/matcher/dfa-classifier.h"
#include "sequence/matcher/seq-matcher.h"
#include "trim-info-simplex.h"

namespace xoos::demux {

namespace sid_scoring {

/**
 * @class SidScoring
 * @brief Shared SID scoring utilities for search-space penalty calculation.
 *
 * Provides Karlin-Altschul-style search-space correction for LUT-based SID matching.
 * Used by independent simplex, paired simplex, and simplex-10x adapters.
 */
class SidScoring {
 public:
  /**
   * @brief Precomputed coefficients for O(1) search-space penalty calculation.
   *
   * For a read of length L, the expected number of random SID hits is represented as
   * `linear * L + constant`. The constant accounts for the shorter number of valid starting
   * positions for longer k-mers near the ends of the read.
   */
  struct SearchSpaceCoefficients {
    /** @brief Coefficient multiplied by read length when estimating random SID hits. */
    f64 linear;
    /** @brief K-mer edge-effect correction added to the length-dependent term. */
    f64 constant;
  };

  /**
   * @brief Computes the per-k probability of a random k-mer matching a SID pool.
   *
   * For each k in [min_k, max_k], calculates unique_kmers / 4^k from the matcher's LUT.
   *
   * @param sid_matcher The SID matcher whose LUT hash tables are queried.
   * @param min_k Minimum k-mer size to evaluate.
   * @param max_k Maximum k-mer size to evaluate.
   * @return Random-match probabilities indexed by `k - min_k`.
   */
  static vec<f64> CalcSIDProbRandomMatch(const SeqMatcher& sid_matcher, u8 min_k, u8 max_k);

  /**
   * @brief Computes the per-k random-match probability across two paired search ends.
   *
   * For each k in [min_k, max_k], calculates `(count_a + count_b) / 4^k` using both SID LUT hash
   * tables. The counts are added rather than deduplicated because paired matching searches both
   * ends and applies one joint correction. This is a conservative upper bound, not a union, and
   * can slightly overestimate the probability when the two pools contain overlapping k-mers.
   * The paired implementation requires matching k-ranges.
   *
   * @param sid_matcher_a SID matcher for the first search end.
   * @param sid_matcher_b SID matcher for the second search end.
   * @param min_k Minimum k-mer size to evaluate.
   * @param max_k Maximum k-mer size to evaluate.
   * @return Combined random-match probabilities indexed by `k - min_k`.
   */
  static vec<f64> CalcSIDProbRandomMatch(const SeqMatcher& sid_matcher_a,
                                         const SeqMatcher& sid_matcher_b,
                                         u8 min_k,
                                         u8 max_k);

  /**
   * @brief Reduces per-k random-match probabilities to linear and constant terms.
   *
   * These coefficients allow CalcSearchSpacePenalty to compute the expected random hit count in
   * O(1), analogous to precomputing the effective search-space parameters in Karlin-Altschul
   * statistics. Expanding `sum_k(prob[k] * (sequence_length - k + 1))` gives:
   *
   * `linear = sum_k(prob[k])`
   *
   * `constant = sum_k(prob[k] * (1 - k))`
   *
   * The expected random hit count is then `linear * sequence_length + constant`.
   *
   * @param sid_kmer_prob Per-k random-match probabilities from CalcSIDProbRandomMatch.
   * @param min_k Minimum k-mer size represented by the first probability.
   * @param max_k Maximum k-mer size represented by the last probability.
   * @return Precomputed linear and constant search-space terms.
   */
  static SearchSpaceCoefficients CalcSearchSpaceCoefficients(const vec<f64>& sid_kmer_prob, u8 min_k, u8 max_k);

  /**
   * @brief Builds search-space coefficients for one independent SID pool.
   *
   * The matcher's SID length and maximum edit distance determine the inclusive k-mer range.
   *
   * @param sid_matcher SID matcher whose LUT defines the independent search space.
   * @return Search-space coefficients for that matcher.
   * @throws error::Error If the matcher's SID length is shorter than its maximum edit distance.
   */
  static SearchSpaceCoefficients MakeSearchSpaceCoefficients(const SeqMatcher& sid_matcher);

  /**
   * @brief Adds independent search-space models for a joint two-end score.
   *
   * Independent 5' and 3' LUT searches contribute distinct random-hit opportunities. Their
   * expected-hit functions are linear in read length, so the joint model is their component-wise
   * sum and supports different SID lengths and edit-distance ranges.
   */
  static SearchSpaceCoefficients AddSearchSpaceCoefficients(const SearchSpaceCoefficients& first,
                                                            const SearchSpaceCoefficients& second);

  /**
   * @brief Builds one joint search-space model for paired 5' and 3' SID observations.
   *
   * Paired observations are scored together, so the random-match opportunities from both ends are
   * combined into one penalty. Both matchers must represent the same inclusive k-mer range.
   *
   * @param sid_matcher_5p SID matcher used at the 5' end.
   * @param sid_matcher_3p SID matcher used at the 3' end.
   * @return Search-space coefficients covering both paired search ends.
   * @throws error::Error If a SID length is invalid or the two matchers have different k-ranges.
   */
  static SearchSpaceCoefficients MakePairedSearchSpaceCoefficients(const SeqMatcher& sid_matcher_5p,
                                                                   const SeqMatcher& sid_matcher_3p);

  /**
   * @brief Computes a length-dependent bit-score penalty for a read.
   *
   * Conceptually equivalent to a Karlin-Altschul E-value correction, this returns
   * `ceil(log2(expected_random_hits))`, floored at zero. Adding the penalty to the base score
   * threshold compensates for the increased number of random matches on longer reads and
   * maintains a consistent false-discovery rate per read.
   *
   * @param coefficients Precomputed coefficients for the relevant SID search space.
   * @param seq_len Length of the read sequence in bases.
   * @return Nonnegative bit-score penalty for the supplied read length.
   */
  static s32 CalcSearchSpacePenalty(const SearchSpaceCoefficients& coefficients, u32 seq_len);

  /**
   * @brief Validates that a SID length is not shorter than its maximum edit distance.
   *
   * @param sid_matcher SID matcher whose configured bounds are validated.
   * @throws error::Error If `length - maximum_edit_distance` would underflow.
   */
  static void ValidateSidBounds(const SeqMatcher& sid_matcher);

  /**
   * @brief Computes the minimum LUT k-mer size for a SID matcher.
   *
   * @param sid_matcher SID matcher whose lower k-mer bound is requested.
   * @return SID length minus maximum edit distance.
   * @throws error::Error If the SID length is shorter than the maximum edit distance.
   */
  static u8 CalcSidMinK(const SeqMatcher& sid_matcher);

  /**
   * @brief Computes the maximum LUT k-mer size for a SID matcher.
   *
   * @param sid_matcher SID matcher whose upper k-mer bound is requested.
   * @return SID length plus maximum edit distance.
   */
  static u8 CalcSidMaxK(const SeqMatcher& sid_matcher);
};

}  // namespace sid_scoring

/**
 * @class DemuxAndTrimSimplex
 * @brief Demultiplexes simplex reads whose 5' and 3' SIDs are statistically independent.
 *
 * Each end is searched, scored, and thresholded against its own SID pool. Only the best candidate
 * per end is retained; there is no cross-end score rescue or per-SID candidate table. When both
 * ends pass independently, their SID pair is mapped to a sample by DualSidMatchingPolicy.
 *
 * This class also owns the seed-extension and scoring primitives shared with paired simplex.
 */
class DemuxAndTrimSimplex {
 public:
  /**
   * @brief Constructs the independent simplex engine.
   *
   * Initializes per-end SID matchers, fixed-region DFA classifiers, anchor searchers, and separate
   * search-space models. The 5' and 3' SID configurations may have different lengths and edit
   * distances because neither end contributes evidence to the other end's acceptance decision.
   *
   * @param param Runtime score and read-length parameters.
   * @param components Adapter components in molecule read order: anchor, SID, stem, stem, SID,
   *   anchor. Components are addressed by position and fixed literals are compiled into the DFA
   *   classifiers and Bitap searchers during construction.
   * @param dual_sid_matching_policy Maps independently accepted SID pairs and permitted partials to samples.
   * @throws error::Error If the component layout or SID bounds are invalid, or if the configured
   *   minimum read length cannot satisfy an independent score threshold.
   */
  DemuxAndTrimSimplex(const DemuxAndTrimParam& param,
                      const vec<AdapterComponent>& components,
                      DualSidMatchingPolicy dual_sid_matching_policy);
  virtual ~DemuxAndTrimSimplex() = default;

  /**
   * @brief Finds the best candidate at each end, applies independent thresholds, and trims the read.
   *
   * The best 5' and 3' candidates are selected without sharing scores or candidate state. Each end
   * must clear its own search-space-adjusted threshold before it can participate in sample mapping.
   * Accepted trim coordinates and SID assignments are written to `record.trim_info_simplex`.
   *
   * @param record Read record to demultiplex and trim; its simplex trim information is updated.
   * @return The TrimInfoSimplex stored in @p record.
   */
  virtual const TrimInfoSimplex& operator()(FixedReadRecord& record) const;

 protected:
  /**
   * @brief Initializes shared simplex matching state with a caller-selected stem tolerance.
   *
   * Derived adapter types use this overload when their fixed stem lengths require a different DFA
   * extension band while retaining the same seed extension and scoring primitives.
   *
   * @param param Runtime score and read-length parameters.
   * @param components Adapter components in molecule order: anchor, SID, stem, stem, SID, anchor.
   * @param dual_sid_matching_policy Policy for mapping independently accepted SID observations.
   * @param stem_max_drop Maximum score drop retained during both insert-facing stem extensions.
   * @throws error::Error If the component layout, SID bounds, or score configuration is invalid.
   */
  DemuxAndTrimSimplex(const DemuxAndTrimParam& param,
                      const vec<AdapterComponent>& components,
                      DualSidMatchingPolicy dual_sid_matching_policy,
                      s32 stem_max_drop);

  /**
   * @brief Backward and forward fixed-region extension results for one SID seed.
   *
   * The two results cover the fixed adapter sequence on either side of a candidate SID match and
   * are combined with the seed's own likelihood score.
   */
  struct PairedDFAResults {
    /** @brief Classification of sequence preceding the SID seed. */
    DFAClassifier::DFAResults bw_dfa_results;
    /** @brief Classification of sequence following the SID seed. */
    DFAClassifier::DFAResults fw_dfa_results;
  };

  /**
   * @brief A SID seed after fixed-region extension and scoring.
   *
   * This compact value is sufficient for candidate comparison and final trim reconstruction in
   * both independent and paired simplex implementations.
   */
  struct ScoredMatch {
    /** @brief Original SID seed match, including barcode identity and edit distance. */
    MatchInfo match;
    /** @brief Log-likelihood score after incorporating both fixed-region extensions. */
    s32 score;
    /** @brief Insert-facing trim boundary derived from the inner stem extension. */
    s32 trim_pos;
  };

  /**
   * @brief Per-read search limits and independent score thresholds.
   *
   * The thresholds are computed separately because the two SID pools can have different sequence
   * lengths, edit distances, LUT contents, and therefore different random-match search spaces.
   */
  struct SearchBounds {
    /** @brief Read length in bases. */
    s32 length;
    /** @brief Minimum number of bases that must remain between trim boundaries to keep the insert non-empty. */
    s32 min_insert_len;
    /** @brief Minimum accepted 5' score: base floor plus the 5' pool search-space penalty. */
    s32 min_score_5p;
    /** @brief Minimum accepted 3' score: base floor plus the 3' pool search-space penalty. */
    s32 min_score_3p;
  };

  /** @brief Runtime demultiplexing and trimming parameters shared by all read evaluations. */
  const DemuxAndTrimParam& _param;

  // The variable SID components are retained because they are referenced for every read. Fixed
  // literals are consumed during construction and baked into the classifiers and searchers below.
  /** @brief SID matcher used to generate candidates from the 5' search region. */
  const SeqMatcherPtr _sid_fwd_matcher;
  /** @brief SID matcher used to generate candidates from the 3' search region. */
  const SeqMatcherPtr _sid_rev_matcher;
  /**
   * @brief Maps independently accepted per-end SID IDs to samples and controls partial assignment.
   */
  const DualSidMatchingPolicy _dual_sid_matching_policy;

  /** @brief Nominal 5' SID length used by MatchInfo::LogLikelihoodScore. */
  const s32 _nominal_seed_size_5p;
  /** @brief Nominal 3' SID length used by MatchInfo::LogLikelihoodScore. */
  const s32 _nominal_seed_size_3p;

  // Fixed literals are compiled into these classifiers at construction and are not otherwise retained.
  /** @brief Backward extension through the outer fixed region left of the 5' SID. */
  const DFAClassifier _bw_extend_anchor_fwd_dfa;
  /** @brief Forward extension through the inner fixed stem right of the 5' SID. */
  const DFAClassifier _fw_extend_stem_fwd_dfa;
  /** @brief Backward extension through the inner fixed stem left of the 3' SID. */
  const DFAClassifier _bw_extend_stem_rev_dfa;
  /** @brief Forward extension through the outer fixed region right of the 3' SID. */
  const DFAClassifier _fw_extend_anchor_rev_dfa;

  /** @brief Forward Bitap searcher for baiting complete searches from the 5' anchor. */
  const Bitap<4> _fw_search_anchor_fwd_bitap;
  /** @brief Backward Bitap searcher for baiting complete searches from the 3' anchor. */
  const Bitap<4> _bw_search_anchor_rev_bitap;

  /** @brief Full-stem matchers used to complete one-SID dual-index partial reads. */
  const Bitap<4> _fw_partial_stem_fwd_bitap;
  const Bitap<4> _bw_partial_stem_fwd_bitap;
  const Bitap<4> _fw_partial_stem_rev_bitap;
  const Bitap<4> _bw_partial_stem_rev_bitap;

  /**
   * @brief Search-space coefficients built solely from the 5' SID pool and k-mer range.
   *
   * These coefficients produce the threshold used only for accepting a 5' observation in the
   * independent algorithm.
   */
  const sid_scoring::SidScoring::SearchSpaceCoefficients _sid_5p_search_space_coefficients;
  /**
   * @brief Search-space coefficients built solely from the 3' SID pool and k-mer range.
   *
   * These coefficients produce the threshold used only for accepting a 3' observation in the
   * independent algorithm.
   */
  const sid_scoring::SidScoring::SearchSpaceCoefficients _sid_3p_search_space_coefficients;

  /**
   * @brief Builds search limits and independent per-end thresholds for one read.
   *
   * @param length Read length in bases.
   * @return Search bounds containing the insert constraint and both adjusted score thresholds.
   */
  SearchBounds MakeSearchBounds(s32 length) const;

  /**
   * @brief Scores one 5' SID seed with its adjacent fixed-region extensions.
   *
   * @param record Read containing the candidate.
   * @param match SID seed match to extend and score.
   * @param bounds Per-read geometry and independent thresholds.
   * @return A scored candidate with its insert start, or std::nullopt if extension geometry is invalid.
   */
  std::optional<ScoredMatch> Score5pMatch(const FixedReadRecord& record,
                                          const MatchInfo& match,
                                          const SearchBounds& bounds) const;

  /**
   * @brief Scores one 3' SID seed with its adjacent fixed-region extensions.
   *
   * @param record Read containing the candidate.
   * @param match SID seed match to extend and score.
   * @param bounds Per-read geometry and independent thresholds.
   * @return A scored candidate with its insert end, or std::nullopt if extension geometry is invalid.
   */
  std::optional<ScoredMatch> Score3pMatch(const FixedReadRecord& record,
                                          const MatchInfo& match,
                                          const SearchBounds& bounds) const;

  /**
   * @brief Searches the 5' region and retains only its highest-scoring valid candidate.
   *
   * @param record Read to search.
   * @param bounds Per-read geometry and independent thresholds.
   * @param search_limit Exclusive upper coordinate for SID seed discovery.
   * @return The best scored 5' candidate, or std::nullopt when none is valid.
   */
  std::optional<ScoredMatch> FindBest5p(const FixedReadRecord& record,
                                        const SearchBounds& bounds,
                                        s32 search_limit) const;

  /**
   * @brief Searches the 3' region and retains only its highest-scoring valid candidate.
   *
   * @param record Read to search.
   * @param bounds Per-read geometry and independent thresholds.
   * @param search_limit Exclusive lower coordinate for SID seed discovery.
   * @return The best scored 3' candidate, or std::nullopt when none is valid.
   */
  std::optional<ScoredMatch> FindBest3p(const FixedReadRecord& record,
                                        const SearchBounds& bounds,
                                        s32 search_limit) const;

  /** @brief Return the inclusive coordinates of a full 5' or 3' stem matched near its read end. */
  std::optional<LociRange> FindTerminalStem(const FixedReadRecord& record, bool is_5p) const;

  /** @brief Return an end-anchored 16bp-or-longer prefix or suffix of a cutoff missing-side stem. */
  std::optional<LociRange> FindTerminalPartialStem(const FixedReadRecord& record, bool is_5p) const;

  /** @brief Complete the missing trim boundary of a uniquely assigned one-SID dual-index read. */
  void CompletePartialStem(const FixedReadRecord& record, const SearchBounds& bounds, TrimInfoSimplex& result) const;

  /**
   * @brief Extends a SID seed in both directions through its adjacent fixed regions.
   *
   * @param record Read containing the SID seed.
   * @param match_info SID seed match that defines the extension boundaries.
   * @param bw_dfa Classifier for sequence preceding the seed.
   * @param fw_dfa Classifier for sequence following the seed.
   * @return Classification results for both fixed-region extensions.
   */
  static PairedDFAResults ExtendCandidate(const FixedReadRecord& record,
                                          const MatchInfo& match_info,
                                          const DFAClassifier& bw_dfa,
                                          const DFAClassifier& fw_dfa);

  /**
   * @brief Validates that the minimum read length can reach the configured score threshold.
   *
   * @param param Runtime parameters containing the base score and minimum read length.
   * @param coefficients Search-space coefficients used to adjust the threshold.
   * @param scoring_mode Human-readable scoring mode included in validation errors.
   * @throws error::Error If the adjusted threshold is unreachable at the minimum read length.
   */
  static void ValidateMinReadLength(const DemuxAndTrimParam& param,
                                    const sid_scoring::SidScoring::SearchSpaceCoefficients& coefficients,
                                    std::string_view scoring_mode);

  /**
   * @brief Maps an independently accepted 5'/3' SID pair to a sample.
   *
   * @param sid_5p_id Barcode ID accepted at the 5' end.
   * @param sid_3p_id Barcode ID accepted at the 3' end.
   * @return The mapped sample ID, or std::nullopt when the pair is not assigned.
   */
  std::optional<u32> PairedSampleId(u32 sid_5p_id, u32 sid_3p_id) const;

  /**
   * @brief Maps a 5'-only match when that partial assignment is unique and permitted.
   *
   * @param sid_5p_id Barcode ID accepted at the 5' end.
   * @return The mapped sample ID, or std::nullopt when partial assignment is unavailable.
   */
  std::optional<u32> SampleIdFor5pPartialMatch(u32 sid_5p_id) const;

  /**
   * @brief Maps a 3'-only match when that partial assignment is unique and permitted.
   *
   * @param sid_3p_id Barcode ID accepted at the 3' end.
   * @return The mapped sample ID, or std::nullopt when partial assignment is unavailable.
   */
  std::optional<u32> SampleIdFor3pPartialMatch(u32 sid_3p_id) const;

  /**
   * @brief Applies independent thresholds and reconstructs the final trim information.
   *
   * Both candidates are filtered using only their own end's threshold. Accepted observations are
   * then passed to the matching policy for full-pair or permitted partial sample assignment.
   *
   * @param bounds Per-read geometry and independent thresholds.
   * @param best_5p Highest-scoring valid 5' candidate, if one was found.
   * @param best_3p Highest-scoring valid 3' candidate, if one was found.
   * @param result Trim information populated with accepted SID and insert coordinates.
   */
  void ReconstructIndependentResult(const SearchBounds& bounds,
                                    const std::optional<ScoredMatch>& best_5p,
                                    const std::optional<ScoredMatch>& best_3p,
                                    TrimInfoSimplex& result) const;
};

}  // namespace xoos::demux
