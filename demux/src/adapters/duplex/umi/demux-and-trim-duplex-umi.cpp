#include "demux-and-trim-duplex-umi.h"

#include "io/read-record.h"

namespace xoos::demux {

// Helper function to create a TrimDuplexUmi for UMI processing
static TrimDuplexUmi CreateTrimUmi(const LutBundleDuplexUmi& lut_bundle) {
  const auto& loop_sequence = lut_bundle.LoopSequence();
  return TrimDuplexUmi{lut_bundle.Sid5pMatcher(),        lut_bundle.Sid3pMatcher(),  lut_bundle.StartSequenceMatcher(),
                       lut_bundle.StopSequenceMatcher(), lut_bundle.Stem5pMatcher(), lut_bundle.Stem3pMatcher(),
                       lut_bundle.Umi5pMatcher(),        lut_bundle.Umi3pMatcher(),  loop_sequence};
}

DemuxAndTrimDuplexUmi::DemuxAndTrimDuplexUmi(const LutBundleDuplexUmi& lut_bundle)
    : DemuxAndTrimDuplexStem(lut_bundle), _trim_umi(CreateTrimUmi(lut_bundle)) {}

std::pair<s32, s32> DemuxAndTrimDuplexUmi::FindUMIPos(FixedReadRecord& read) const {
  // Process the read and UMI using the TrimDuplexUmi instance
  return _trim_umi.FindUMI(read);
}

}  // namespace xoos::demux
