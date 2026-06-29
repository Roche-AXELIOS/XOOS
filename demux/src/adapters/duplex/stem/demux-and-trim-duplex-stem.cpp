#include "demux-and-trim-duplex-stem.h"

#include "io/read-record.h"

namespace xoos::demux {

// Helper function to create a TrimDuplexStem for stem processing
TrimDuplexStem CreateTrimStem(const LutBundleDuplexStem& lut_bundle) {
  return TrimDuplexStem{lut_bundle.Sid5pMatcher(),        lut_bundle.Sid3pMatcher(),  lut_bundle.StartSequenceMatcher(),
                        lut_bundle.StopSequenceMatcher(), lut_bundle.Stem5pMatcher(), lut_bundle.Stem3pMatcher(),
                        lut_bundle.LoopSequence()};
}

DemuxAndTrimDuplexStem::DemuxAndTrimDuplexStem(const LutBundleDuplexStem& lut_bundle)
    : DemuxAndTrimDuplex(lut_bundle), _trim_stem(CreateTrimStem(lut_bundle)) {}

std::pair<s32, s32> DemuxAndTrimDuplexStem::FindUMIPos(FixedReadRecord& read) const {
  // Process the read and find stem positions using the TrimDuplexStem instance
  return _trim_stem.FindStem(read);
}

}  // namespace xoos::demux
