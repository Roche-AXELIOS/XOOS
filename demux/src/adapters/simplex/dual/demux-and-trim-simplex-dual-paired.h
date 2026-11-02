#pragma once

#include <limits>
#include <optional>

#include <gtl/phmap.hpp>

#include "adapters/simplex/demux-and-trim-simplex.h"

namespace xoos::demux {

/// @brief Demultiplexes dual-index simplex reads using paired sample-sheet evidence.
///
/// Unlike normal paired simplex, the SID identities at the two ends are different. The sample
/// sheet defines the valid pairs, so this class retains the best 5' candidate for each observed
/// SID and joins a 3' candidate only to its legal 5' partners. The resulting two-end score is
/// accepted against one joint search-space threshold. When no legal pair is accepted, the
/// existing independent reconstruction preserves dual-sid-mode partial-read behavior.
class DemuxAndTrimSimplexDualPaired final : public DemuxAndTrimSimplex {
 public:
  DemuxAndTrimSimplexDualPaired(const DemuxAndTrimParam& param,
                                const vec<AdapterComponent>& components,
                                DualSidMatchingPolicy dual_sid_matching_policy);

  const TrimInfoSimplex& operator()(FixedReadRecord& record) const override;

 private:
  using Candidates5p = gtl::flat_hash_map<u32, ScoredMatch>;

  struct Search5pResult {
    Candidates5p candidates;
    std::optional<ScoredMatch> best;
  };

  struct Search3pResult {
    std::optional<ScoredMatch> best;
    s32 best_pair_score{std::numeric_limits<s32>::min()};
  };

  const sid_scoring::SidScoring::SearchSpaceCoefficients _paired_search_space_coefficients;

  Search5pResult FindPaired5pCandidates(const FixedReadRecord& record,
                                        const SearchBounds& bounds,
                                        s32 search_limit) const;

  Search3pResult FindPaired3pCandidates(const FixedReadRecord& record,
                                        const SearchBounds& bounds,
                                        s32 search_limit,
                                        s32 min_pair_score,
                                        const Candidates5p& candidates_5p,
                                        TrimInfoSimplex& result) const;
};

}  // namespace xoos::demux
