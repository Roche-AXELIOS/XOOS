#include "demux-and-trim-simplex-dual-paired.h"

#include <utility>

#include <xoos/types/vec.h>

namespace xoos::demux {

DemuxAndTrimSimplexDualPaired::DemuxAndTrimSimplexDualPaired(const DemuxAndTrimParam& param,
                                                             const vec<AdapterComponent>& components,
                                                             DualSidMatchingPolicy dual_sid_matching_policy)
    : DemuxAndTrimSimplex{param, components, std::move(dual_sid_matching_policy)},
      _paired_search_space_coefficients{sid_scoring::SidScoring::AddSearchSpaceCoefficients(
          _sid_5p_search_space_coefficients, _sid_3p_search_space_coefficients)} {
  ValidateMinReadLength(param, _paired_search_space_coefficients, "dual paired SID");
}

DemuxAndTrimSimplexDualPaired::Search5pResult DemuxAndTrimSimplexDualPaired::FindPaired5pCandidates(
    const FixedReadRecord& record, const SearchBounds& bounds, const s32 search_limit) const {
  Search5pResult result;
  auto search_pos = 0;
  vec<MatchInfo> matches;
  while (search_pos < search_limit) {
    _sid_fwd_matcher->FindNextBarcode(
        ReadEnd::k5p, search_pos, record.TwoBitsSeq(), static_cast<size_t>(search_limit), matches);
    if (matches.empty()) {
      break;
    }
    search_pos = static_cast<s32>(matches.front().SPos()) + 1;
    for (const auto& match : matches) {
      const auto candidate = Score5pMatch(record, match, bounds);
      if (!candidate.has_value()) {
        continue;
      }
      const auto sid_5p_id = candidate->match.BarcodeId();
      const auto existing = result.candidates.find(sid_5p_id);
      if (existing == result.candidates.end() || candidate->score > existing->second.score) {
        result.candidates.insert_or_assign(sid_5p_id, *candidate);
      }
      if (!result.best.has_value() || candidate->score > result.best->score) {
        result.best = candidate;
      }
    }
  }
  return result;
}

DemuxAndTrimSimplexDualPaired::Search3pResult DemuxAndTrimSimplexDualPaired::FindPaired3pCandidates(
    const FixedReadRecord& record,
    const SearchBounds& bounds,
    const s32 search_limit,
    const s32 min_pair_score,
    const Candidates5p& candidates_5p,
    TrimInfoSimplex& result) const {
  Search3pResult search_result;
  auto search_pos = bounds.length;
  vec<MatchInfo> matches;
  while (search_pos > search_limit) {
    _sid_rev_matcher->FindNextBarcode(
        ReadEnd::k3p, search_pos, record.TwoBitsSeq(), static_cast<size_t>(bounds.length), matches);
    if (matches.empty()) {
      break;
    }
    search_pos = static_cast<s32>(matches.front().EPos()) - 1;
    for (const auto& match : matches) {
      const auto candidate_3p = Score3pMatch(record, match, bounds);
      if (!candidate_3p.has_value()) {
        continue;
      }
      if (!search_result.best.has_value() || candidate_3p->score > search_result.best->score) {
        search_result.best = candidate_3p;
      }

      for (const auto& pairing : _dual_sid_matching_policy.PairingsFor3pSid(match.BarcodeId())) {
        const auto candidate_5p = candidates_5p.find(pairing.sid_5p_id);
        if (candidate_5p == candidates_5p.end()) {
          continue;
        }
        const auto pair_score = candidate_5p->second.score + candidate_3p->score;
        if (pair_score < min_pair_score || pair_score <= search_result.best_pair_score ||
            candidate_5p->second.trim_pos >= candidate_3p->trim_pos) {
          continue;
        }

        result.sid = pairing.sample_id;
        result.sid_5p = pairing.sample_id;
        result.sid_5p_edist = candidate_5p->second.match.EDist();
        result.score_5p = candidate_5p->second.score;
        result.sid_3p = pairing.sample_id;
        result.sid_3p_edist = candidate_3p->match.EDist();
        result.score_3p = candidate_3p->score;
        result.insert.spos = static_cast<u32>(candidate_5p->second.trim_pos);
        result.insert.epos = static_cast<u32>(candidate_3p->trim_pos);
        search_result.best_pair_score = pair_score;
      }
    }
  }
  return search_result;
}

const TrimInfoSimplex& DemuxAndTrimSimplexDualPaired::operator()(FixedReadRecord& record) const {
  auto& result = record.trim_info_simplex;
  result.Clear();
  const auto bounds = MakeSearchBounds(static_cast<s32>(record.SeqLen()));
  result.insert.epos = static_cast<u32>(bounds.length);
  const auto search_5p = FindPaired5pCandidates(record, bounds, bounds.length - bounds.min_insert_len);
  const auto min_pair_score =
      _param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(_paired_search_space_coefficients,
                                                                         static_cast<u32>(bounds.length));
  const auto search_3p =
      FindPaired3pCandidates(record, bounds, bounds.min_insert_len, min_pair_score, search_5p.candidates, result);
  if (search_3p.best_pair_score >= min_pair_score) {
    return result;
  }

  result.Clear();
  ReconstructIndependentResult(bounds, search_5p.best, search_3p.best, result);
  CompletePartialStem(record, bounds, result);
  return result;
}

}  // namespace xoos::demux
