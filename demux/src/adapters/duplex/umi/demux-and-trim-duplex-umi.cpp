#include "demux-and-trim-duplex-umi.h"

#include "io/read-record.h"

namespace xoos::demux {

// Helper function to create a TrimDuplexUmi for UMI processing
static TrimDuplexUmi CreateTrimUmi(const LutBundleDuplexUmi& lut_bundle) {
  return TrimDuplexUmi{lut_bundle.Runway1(),
                       lut_bundle.StartSequence(),
                       lut_bundle.Stem4(),
                       lut_bundle.Umi6Matcher(),
                       lut_bundle.Umi10Matcher()};
}

DemuxAndTrimDuplexUmi::DemuxAndTrimDuplexUmi(const LutBundleDuplexUmi& lut_bundle)
    : DemuxAndTrimDuplexStem(lut_bundle), _trim_umi(CreateTrimUmi(lut_bundle)) {
}

std::pair<s32, s32> DemuxAndTrimDuplexUmi::FindUMIPos(FixedReadRecord& read) const {
  // Process the read and UMI using the TrimDuplexUmi instance
  return _trim_umi.FindUMI(read);
}

}  // namespace xoos::demux
