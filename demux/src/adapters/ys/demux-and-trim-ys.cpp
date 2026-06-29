#include "demux-and-trim-ys.h"

#include <algorithm>

#include "io/read-record.h"

namespace xoos::demux {
static Trim5pYs CreateTrim5p(const bool enable_partial, const LutBundleYs& lut_bundle) {
  return Trim5pYs{enable_partial, lut_bundle.Runway5pMatcher(), lut_bundle.Sid5pMatcher(),
                  lut_bundle.SidSpacer5pMatcher()};
}

static Trim3pYs CreateTrim3p(const bool enable_partial, const LutBundleYs& lut_bundle) {
  return Trim3pYs{enable_partial, lut_bundle.Sid3pMatcher(), lut_bundle.SidSpacer3pMatcher()};
}

DemuxAndTrimYs::DemuxAndTrimYs(const bool enable_partial, const LutBundleYs& lut_bundle)
    : _enable_partial(enable_partial),
      _trim_5p{CreateTrim5p(enable_partial, lut_bundle)},
      _trim_3p{CreateTrim3p(enable_partial, lut_bundle)} {}

TrimInfoYs DemuxAndTrimYs::operator()(const FixedReadRecord& record) const {
  // speed optimization by converting the sequence to a 2-bit representation.
  const auto trim_5p = _trim_5p.Trim(record.TwoBitsSeq(), record.SeqLen(), record.Seq());
  const auto trim_3p = _trim_3p.Trim(record.TwoBitsSeq(), record.SeqLen(), record.Seq(), trim_5p.insert_start);
  return Demux(trim_5p, trim_3p);
}

TrimInfoYs DemuxAndTrimYs::Demux(const Trim5pInfoYs& trim_5p, const Trim3pInfoYs& trim_3p) const {
  const auto insert_start = trim_5p.insert_start;
  const auto insert_end = trim_3p.insert_end;
  return TrimInfoYs{
      DetermineSampleId(trim_5p, trim_3p), MatchInfo::BarcodeId(trim_5p.sid_match),
      MatchInfo::EDist(trim_5p.sid_match), MatchInfo::BarcodeId(trim_3p.sid_match),
      MatchInfo::EDist(trim_3p.sid_match), LociRange{std::min(insert_start, insert_end), insert_end},
  };
}

std::optional<u32> DemuxAndTrimYs::DetermineSampleId(const Trim5pInfoYs& trim_5p, const Trim3pInfoYs& trim_3p) const {
  const auto sid_edist_5p = MatchInfo::EDist(trim_5p.sid_match);
  const auto sid_edist_3p = MatchInfo::EDist(trim_3p.sid_match);
  if (!_enable_partial && !(trim_5p.sid_match && trim_3p.sid_match)) {
    return std::nullopt;
  }

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
