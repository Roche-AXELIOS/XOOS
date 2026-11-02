#pragma once

#include <limits>
#include <optional>

#include <gtl/phmap.hpp>

#include "adapters/simplex/demux-and-trim-simplex.h"

namespace xoos::demux {

/**
 * @class DemuxAndTrimSimplexPaired
 * @brief Demultiplexes normal YS reads whose 5' and 3' observations represent the same SID.
 *
 * Identifies barcode SIDs at both ends using LUT seed matching followed by DFA-based extension
 * through the adjacent fixed adapter regions. Each candidate receives an alignment-like
 * log-likelihood score, and concordant observations of the same SID can contribute combined
 * evidence to the final decision and insert boundaries.
 *
 * Paired matching retains the best 5' candidate for every SID so a later 3' observation can form
 * a concordant pair even when that SID was not the globally best 5' candidate. A pair is evaluated
 * with combined 5'/3' evidence against one joint random-match threshold. The two matchers must have
 * the same SID length and maximum edit distance.
 */
class DemuxAndTrimSimplexPaired final : public DemuxAndTrimSimplex {
 public:
  /**
   * @brief Constructs the paired simplex demultiplexing and trimming engine.
   *
   * Initializes the shared SID matchers, fixed-region DFA classifiers, anchor Bitap searchers, and
   * the joint random-match correction used for paired scoring. The SID pools must use compatible
   * lengths and edit-distance ranges because candidates are joined by one barcode identity and
   * evaluated against one combined threshold.
   *
   * @param param Runtime demultiplexing and trimming parameters, including score and read-length
   *   limits.
   * @param components Adapter components in molecule read order: anchor, SID, stem, stem, SID,
   *   anchor. Components are addressed by position in the shared simplex constructor.
   * @throws error::Error If the component layout or SID bounds are invalid, the two SID matchers
   *   have incompatible k-mer ranges, or the minimum read length cannot reach the joint threshold.
   */
  DemuxAndTrimSimplexPaired(const DemuxAndTrimParam& param, const vec<AdapterComponent>& components);

  /**
   * @brief Finds the best concordant pair or falls back to the strongest accepted single end.
   *
   * Searches for 5' and 3' SID seeds, extends each seed through its flanking fixed regions, and
   * scores the resulting candidates by log likelihood. Same-SID candidates are tested using their
   * combined score. If no concordant pair outranks the singleton evidence, the strongest accepted
   * single end determines the SID and the available trim boundary. Results are written to
   * `record.trim_info_simplex`.
   *
   * @param record Read record to demultiplex and trim; its simplex trim information is updated.
   * @return A const reference to the updated TrimInfoSimplex stored in @p record.
   */
  const TrimInfoSimplex& operator()(FixedReadRecord& record) const override;

 private:
  /**
   * @brief Best 5' candidate for each SID ID.
   *
   * Keeping one candidate per barcode is required because the globally best 5' candidate may not
   * be the SID observed at 3'. The table allows every valid 3' candidate to look up its concordant
   * 5' evidence directly.
   */
  using Candidates5p = gtl::flat_hash_map<u32, ScoredMatch>;

  /**
   * @brief Output of the 5' search stage.
   *
   * Contains the per-SID table needed by the pairing stage and the globally strongest 5' candidate
   * needed to compare paired evidence against the best singleton.
   */
  struct PairedSearch5pResult {
    /** @brief Highest-scoring valid 5' candidate for every SID ID. */
    Candidates5p candidates;
    /** @brief Highest-scoring valid 5' candidate across all SID IDs. */
    std::optional<ScoredMatch> best;
  };

  /**
   * @brief Output of the 3' search and pairing stage.
   *
   * The best unpaired 3' observation is retained for singleton reconstruction. The best pair score
   * tracks whether the TrimInfoSimplex populated during pairing outranks both singleton choices.
   */
  struct PairedSearch3pResult {
    /** @brief Highest-scoring 3' candidate that had no concordant 5' SID candidate. */
    std::optional<ScoredMatch> best_unpaired;
    /** @brief Highest combined score among accepted concordant pairs. */
    s32 best_pair_score{std::numeric_limits<s32>::min()};
  };

  /**
   * @brief Joint search-space model covering the paired 5' and 3' search opportunities.
   *
   * The expected random hit count has the form `linear * read_length + constant`, with both SID
   * LUTs contributing to the coefficients. Its length-dependent penalty is added once to the base
   * threshold used for combined pair scores and paired-mode singleton fallback.
   */
  const sid_scoring::SidScoring::SearchSpaceCoefficients _paired_search_space_coefficients;

  /**
   * @brief Searches the 5' side and collects the candidates needed for later pairing.
   *
   * Every SID keeps its highest-scoring valid candidate, while a second value tracks the strongest
   * 5' candidate across the entire pool for possible singleton reconstruction.
   *
   * @param record Read being processed.
   * @param bounds Read-level geometry and score thresholds prepared by the base class.
   * @param search_limit Exclusive upper coordinate for 5' SID seed discovery.
   * @return Per-SID 5' candidates and the globally strongest 5' candidate.
   */
  PairedSearch5pResult FindPaired5pCandidates(const FixedReadRecord& record,
                                              const SearchBounds& bounds,
                                              s32 search_limit) const;

  /**
   * @brief Searches the 3' region and pairs each candidate with the same SID retained at 5'.
   *
   * Each valid 3' candidate is looked up by barcode ID in the 5' candidate table. A concordant pair
   * must clear the joint threshold and leave an ordered, nonempty insert. The best accepted pair is
   * written directly into @p result, while the strongest 3' candidate without a 5' partner is
   * retained for singleton comparison.
   *
   * @param record Read being processed.
   * @param bounds Read-level geometry and score thresholds prepared by the base class.
   * @param search_limit Exclusive lower coordinate for 3' SID seed discovery.
   * @param min_pair_score Minimum combined score required to accept a concordant pair.
   * @param candidates_5p Per-SID candidates produced by FindPaired5pCandidates.
   * @param result Output trim information updated whenever a stronger accepted pair is found.
   * @return The strongest unpaired 3' candidate and the best accepted pair score.
   */
  PairedSearch3pResult FindPaired3pCandidate(const FixedReadRecord& record,
                                             const SearchBounds& bounds,
                                             s32 search_limit,
                                             s32 min_pair_score,
                                             const Candidates5p& candidates_5p,
                                             TrimInfoSimplex& result) const;

  /**
   * @brief Reconstructs a partial read when no concordant pair wins.
   *
   * Both singleton candidates are tested against the paired-mode threshold. Their per-end SID,
   * edit-distance, and score fields remain available in the result, and the stronger accepted end
   * supplies the sample SID and its corresponding insert boundary.
   *
   * @param bounds Read-level geometry used to compute the adjusted singleton threshold.
   * @param best_5p Strongest valid 5' candidate found during the first search stage.
   * @param best_3p Strongest valid 3' candidate without a concordant 5' partner.
   * @param result Output trim information populated with the accepted singleton evidence.
   */
  void ReconstructPairedPartial(const SearchBounds& bounds,
                                const std::optional<ScoredMatch>& best_5p,
                                const std::optional<ScoredMatch>& best_3p,
                                TrimInfoSimplex& result) const;
};

}  // namespace xoos::demux
