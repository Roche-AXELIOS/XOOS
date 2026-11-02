#include "demux-and-trim-simplex-paired.h"

#include <algorithm>
#include <limits>

#include <xoos/types/vec.h>

namespace xoos::demux {

namespace {
// Short paired stems can use a wider extension band without exceeding the DFA state-space limit.
constexpr s32 kPairedStemMaxDrop = 10;
}  // namespace

DemuxAndTrimSimplexPaired::DemuxAndTrimSimplexPaired(const DemuxAndTrimParam& param,
                                                     const vec<AdapterComponent>& components)
    : DemuxAndTrimSimplex{param, components, DualSidMatchingPolicy{}, kPairedStemMaxDrop},
      _paired_search_space_coefficients{
          sid_scoring::SidScoring::MakePairedSearchSpaceCoefficients(*_sid_fwd_matcher, *_sid_rev_matcher)} {
  ValidateMinReadLength(param, _paired_search_space_coefficients, "paired SID");
}

DemuxAndTrimSimplexPaired::PairedSearch5pResult DemuxAndTrimSimplexPaired::FindPaired5pCandidates(
    const FixedReadRecord& record, const SearchBounds& bounds, const s32 search_limit) const {
  PairedSearch5pResult result;
  result.candidates.reserve(_sid_fwd_matcher->Pool().size());
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
      auto candidate = Score5pMatch(record, match, bounds);
      if (!candidate) {
        continue;
      }
      const auto sid = candidate->match.BarcodeId();
      const auto it = result.candidates.find(sid);
      if (it == result.candidates.end() || candidate->score > it->second.score) {
        result.candidates.insert_or_assign(sid, *candidate);
      }
      if (!result.best || candidate->score > result.best->score) {
        result.best = candidate;
      }
    }
  }
  return result;
}

DemuxAndTrimSimplexPaired::PairedSearch3pResult DemuxAndTrimSimplexPaired::FindPaired3pCandidate(
    const FixedReadRecord& record,
    const SearchBounds& bounds,
    const s32 search_limit,
    const s32 min_pair_score,
    const Candidates5p& candidates_5p,
    TrimInfoSimplex& result) const {
  PairedSearch3pResult search_result;
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
      auto candidate_3p = Score3pMatch(record, match, bounds);
      if (!candidate_3p) {
        continue;
      }
      const auto candidate_5p = candidates_5p.find(match.BarcodeId());
      if (candidate_5p == candidates_5p.end()) {
        if (!search_result.best_unpaired || candidate_3p->score > search_result.best_unpaired->score) {
          search_result.best_unpaired = candidate_3p;
        }
        continue;
      }
      const auto pair_score = candidate_5p->second.score + candidate_3p->score;
      if (pair_score > search_result.best_pair_score && pair_score >= min_pair_score &&
          candidate_5p->second.trim_pos < candidate_3p->trim_pos) {
        result.sid = match.BarcodeId();
        result.sid_5p = match.BarcodeId();
        result.sid_5p_edist = candidate_5p->second.match.EDist();
        result.score_5p = candidate_5p->second.score;
        result.sid_3p = match.BarcodeId();
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

void DemuxAndTrimSimplexPaired::ReconstructPairedPartial(const SearchBounds& bounds,
                                                         const std::optional<ScoredMatch>& best_5p,
                                                         const std::optional<ScoredMatch>& best_3p,
                                                         TrimInfoSimplex& result) const {
  const auto min_score = _param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(
                                                _paired_search_space_coefficients, static_cast<u32>(bounds.length));
  const bool has_5p = best_5p && best_5p->score >= min_score;
  const bool has_3p = best_3p && best_3p->score >= min_score;
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
  if (has_3p && (!has_5p || best_3p->score > best_5p->score)) {
    result.sid = best_3p->match.BarcodeId();
    result.insert.epos = static_cast<u32>(best_3p->trim_pos);
  } else if (has_5p) {
    result.sid = best_5p->match.BarcodeId();
    result.insert.spos = static_cast<u32>(best_5p->trim_pos);
    result.insert.epos = static_cast<u32>(bounds.length);
  }
}

const TrimInfoSimplex& DemuxAndTrimSimplexPaired::operator()(FixedReadRecord& record) const {
  auto& result = record.trim_info_simplex;
  result.Clear();
  const auto bounds = MakeSearchBounds(static_cast<s32>(record.SeqLen()));
  result.insert.epos = static_cast<u32>(bounds.length);
  const auto search_5p = FindPaired5pCandidates(record, bounds, bounds.length - bounds.min_insert_len);
  const auto min_pair_score =
      _param.min_score + sid_scoring::SidScoring::CalcSearchSpacePenalty(_paired_search_space_coefficients,
                                                                         static_cast<u32>(bounds.length));
  const auto search_3p =
      FindPaired3pCandidate(record, bounds, bounds.min_insert_len, min_pair_score, search_5p.candidates, result);
  const auto best_5p_score = search_5p.best ? search_5p.best->score : std::numeric_limits<s32>::min();
  const auto best_3p_score = search_3p.best_unpaired ? search_3p.best_unpaired->score : std::numeric_limits<s32>::min();
  if (search_3p.best_pair_score >= std::max(best_5p_score, best_3p_score)) {
    return result;
  }
  result.Clear();
  ReconstructPairedPartial(bounds, search_5p.best, search_3p.best_unpaired, result);
  return result;
}

}  // namespace xoos::demux
