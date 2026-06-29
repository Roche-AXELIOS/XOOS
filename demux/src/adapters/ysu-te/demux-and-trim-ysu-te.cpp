#include "demux-and-trim-ysu-te.h"

#include <algorithm>

#include "io/read-record.h"

namespace xoos::demux {

Trim5pYsuTe CreateTrim5p(bool enable_partial, const LutBundleYsuTe& lut_bundle) {
  return Trim5pYsuTe{enable_partial, lut_bundle.Runway5pMatcher(), lut_bundle.Sid5pMatcher(),
                     lut_bundle.SidUmiSpacer5pMatcher(), lut_bundle.Umi5pMatcher()};
}

Trim3pYsuTe CreateTrim3p(bool enable_partial, const LutBundleYsuTe& lut_bundle) {
  return Trim3pYsuTe{enable_partial, lut_bundle.Sid3pMatcher(), lut_bundle.SidUmiSpacer3pMatcher(),
                     lut_bundle.Umi3pMatcher()};
}

DemuxAndTrimYsuTe::DemuxAndTrimYsuTe(bool enable_partial, const LutBundleYsuTe& lut_bundle)
    : _enable_partial{enable_partial},
      _trim_5p{CreateTrim5p(enable_partial, lut_bundle)},
      _trim_3p{CreateTrim3p(enable_partial, lut_bundle)} {}

TrimInfoYsuTe DemuxAndTrimYsuTe::operator()(const FixedReadRecord& record) const {
  // speed optimization by converting the sequence to a 2-bit representation.
  auto trim_5p = _trim_5p.Trim(record.TwoBitsSeq(), record.SeqLen(), record.Seq());
  auto trim_3p = _trim_3p.Trim(record.TwoBitsSeq(), record.SeqLen(), record.Seq(), trim_5p.insert_start);
  return Demux(trim_5p, trim_3p);
}

TrimInfoYsuTe DemuxAndTrimYsuTe::Demux(const Trim5pInfoYsuTe& trim_5p, const Trim3pInfoYsuTe& trim_3p) const {
  auto insert_start = trim_5p.insert_start;
  auto insert_end = trim_3p.insert_end;

  return TrimInfoYsuTe{
      DetermineSampleId(trim_5p, trim_3p),     MatchInfo::BarcodeId(trim_5p.sid_match),
      MatchInfo::EDist(trim_5p.sid_match),     MatchInfo::BarcodeId(trim_5p.umi_match),
      MatchInfo::BarcodeId(trim_3p.sid_match), MatchInfo::EDist(trim_3p.sid_match),
      MatchInfo::BarcodeId(trim_3p.umi_match), LociRange{std::min(insert_start, insert_end), insert_end},
  };
}

std::optional<uint> DemuxAndTrimYsuTe::DetermineSampleId(const Trim5pInfoYsuTe& trim_5p,
                                                         const Trim3pInfoYsuTe& trim_3p) const {
  auto has_umi = _enable_partial ? trim_5p.umi_match || trim_3p.umi_match : trim_5p.umi_match && trim_3p.umi_match;
  if (!has_umi) {
    return std::nullopt;
  }

  auto sid_edist_5p = MatchInfo::EDist(trim_5p.sid_match);
  auto sid_edist_3p = MatchInfo::EDist(trim_3p.sid_match);

  if (sid_edist_5p && !sid_edist_3p) {
    return MatchInfo::BarcodeId(trim_5p.sid_match);
  }

  if (!sid_edist_5p && sid_edist_3p) {
    return MatchInfo::BarcodeId(trim_3p.sid_match);
  }

  if (!sid_edist_5p && !sid_edist_3p) {
    return std::nullopt;
  }

  if (sid_edist_3p < sid_edist_5p) {
    return MatchInfo::BarcodeId(trim_3p.sid_match);
  }

  return MatchInfo::BarcodeId(trim_5p.sid_match);
}

}  // namespace xoos::demux
