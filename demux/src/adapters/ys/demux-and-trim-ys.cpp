#include "demux-and-trim-ys.h"

#include "io/read-record.h"

namespace xoos::demux {
static Trim5pYs CreateTrim5p(const LutBundleYs& lut_bundle) {
  return Trim5pYs{lut_bundle.Runway1(), lut_bundle.Sid2Matcher()};
}

static Trim3pYs CreateTrim3p(const LutBundleYs& lut_bundle) {
  return Trim3pYs{lut_bundle.Sid4Matcher(), lut_bundle.Bait5()};
}

DemuxAndTrimYs::DemuxAndTrimYs(const LutBundleYs& lut_bundle)
    : _trim_5p{CreateTrim5p(lut_bundle)}, _trim_3p{CreateTrim3p(lut_bundle)} {
}

TrimInfoSimplex DemuxAndTrimYs::operator()(const FixedReadRecord& record) const {
  // speed optimization by converting the sequence to a 2-bit representation.
  const auto trim_5p = _trim_5p.Trim(record.TwoBitsSeq(), record.SeqLen(), record.Seq());
  const auto trim_3p = _trim_3p.Trim(record.TwoBitsSeq(), record.SeqLen(), record.Seq(), trim_5p.insert_start);
  return Demux(trim_5p, trim_3p, record.SeqLen());
}

TrimInfoSimplex DemuxAndTrimYs::Demux(const Trim5pInfoYs& trim_5p, const Trim3pInfoYs& trim_3p, const u32 seq_length) {
  const auto sid_5p = MatchInfo::BarcodeId(trim_5p.sid_match);
  const auto sid_3p = MatchInfo::BarcodeId(trim_3p.sid_match);
  const auto sid_5p_edist = MatchInfo::EDist(trim_5p.sid_match);
  const auto sid_3p_edist = MatchInfo::EDist(trim_3p.sid_match);

  auto insert_start = trim_5p.insert_start;
  auto insert_end = trim_3p.insert_end;

  // Winning-side-only trimming for discordant reads
  if (sid_5p && sid_3p && *sid_5p != *sid_3p) {
    if (sid_5p_edist && sid_3p_edist && *sid_3p_edist < *sid_5p_edist) {
      // 3' wins — don't trim 5'
      insert_start = 0;
    } else {
      // 5' wins (including tie) — don't trim 3'
      insert_end = seq_length;
    }
  }

  return TrimInfoSimplex{
      DetermineSampleId(trim_5p, trim_3p),
      sid_5p,
      sid_5p_edist,
      sid_3p,
      sid_3p_edist,
      LociRange{std::min(insert_start, insert_end), insert_end},
      // Setting score here is a hack to allow for smoother breaking upstream
      sid_5p_edist,
      sid_3p_edist,
  };
}

std::optional<u32> DemuxAndTrimYs::DetermineSampleId(const Trim5pInfoYs& trim_5p, const Trim3pInfoYs& trim_3p) {
  // Assignment is SID-only: a read is assigned whenever at least one SID is identified.
  const auto sid_edist_5p = MatchInfo::EDist(trim_5p.sid_match);
  const auto sid_edist_3p = MatchInfo::EDist(trim_3p.sid_match);

  if (sid_edist_5p.has_value() && !sid_edist_3p.has_value()) {
    return MatchInfo::BarcodeId(trim_5p.sid_match);
  }

  if (!sid_edist_5p.has_value() && sid_edist_3p.has_value()) {
    return MatchInfo::BarcodeId(trim_3p.sid_match);
  }

  // At this point sid_edist_5p and sid_edist_3p are either both set or both unset, so testing one suffices.
  if (!sid_edist_5p.has_value()) {
    return std::nullopt;
  }

  if (sid_edist_3p < sid_edist_5p) {
    return MatchInfo::BarcodeId(trim_3p.sid_match);
  }

  return MatchInfo::BarcodeId(trim_5p.sid_match);
}

}  // namespace xoos::demux
