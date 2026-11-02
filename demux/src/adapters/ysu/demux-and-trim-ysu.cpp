#include "demux-and-trim-ysu.h"

#include <xoos/types/int.h>

#include "io/read-record.h"

namespace xoos::demux {

static Trim5pYsu CreateTrim5p(const LutBundleYsu& lut_bundle) {
  return Trim5pYsu{lut_bundle.Sid1Matcher(), lut_bundle.Umi2Matcher(), lut_bundle.Bait1()};
}

static Trim3pYsu CreateTrim3p(const LutBundleYsu& lut_bundle) {
  return Trim3pYsu{lut_bundle.Sid4Matcher(), lut_bundle.Umi3Matcher(), lut_bundle.Bait4()};
}

DemuxAndTrimYsu::DemuxAndTrimYsu(const LutBundleYsu& lut_bundle)
    : _trim_5p{CreateTrim5p(lut_bundle)}, _trim_3p{CreateTrim3p(lut_bundle)} {
}

TrimInfoYsu DemuxAndTrimYsu::operator()(const FixedReadRecord& record) const {
  // speed optimization by converting the sequence to a 2-bit representation.
  const auto trim_5p = _trim_5p.Trim(record.TwoBitsSeq(), record.SeqLen(), record.Seq());
  const auto trim_3p = _trim_3p.Trim(record.TwoBitsSeq(), record.SeqLen(), record.Seq(), trim_5p.insert_start);
  return Demux(trim_5p, trim_3p, record.SeqLen());
}

TrimInfoYsu DemuxAndTrimYsu::Demux(const Trim5pInfoYsu& trim_5p, const Trim3pInfoYsu& trim_3p, const u32 seq_length) {
  const auto sid_5p = MatchInfo::BarcodeId(trim_5p.sid_match);
  const auto sid_3p = MatchInfo::BarcodeId(trim_3p.sid_match);
  const auto sid_5p_edist = MatchInfo::EDist(trim_5p.sid_match);
  const auto sid_3p_edist = MatchInfo::EDist(trim_3p.sid_match);

  auto insert_start = trim_5p.insert_start;
  auto insert_end = trim_3p.insert_end;

  auto umi_5p = MatchInfo::BarcodeId(trim_5p.umi_match);
  auto umi_3p = MatchInfo::BarcodeId(trim_3p.umi_match);

  // Winning-side-only trimming for discordant reads: clear the losing side's trim and UMI
  if (sid_5p && sid_3p && *sid_5p != *sid_3p) {
    if (sid_5p_edist && sid_3p_edist && *sid_3p_edist < *sid_5p_edist) {
      // 3' wins — don't trim 5', clear 5' UMI
      insert_start = 0;
      umi_5p = std::nullopt;
    } else {
      // 5' wins (including tie) — don't trim 3', clear 3' UMI
      insert_end = seq_length;
      umi_3p = std::nullopt;
    }
  }

  return TrimInfoYsu{
      DetermineSampleId(trim_5p, trim_3p),
      sid_5p,
      sid_5p_edist,
      umi_5p,
      sid_3p,
      sid_3p_edist,
      umi_3p,
      LociRange{std::min(insert_start, insert_end), insert_end},
  };
}

std::optional<u32> DemuxAndTrimYsu::DetermineSampleId(const Trim5pInfoYsu& trim_5p, const Trim3pInfoYsu& trim_3p) {
  // Assignment is SID-only: a read is assigned to a sample whenever at least one SID is identified,
  // regardless of UMI presence. The UMI requirement is enforced downstream as a post-assignment
  // filter (the umi_failed gate), so that UMI-less reads are counted as failed_assigned rather than
  // unassigned.
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
