#include "demux-and-trim-duplex.h"

#include "io/read-record.h"

namespace xoos::demux {

struct DuplexMetrics;

static TrimDuplex CreateTrim(const LutBundleDuplex& lut_bundle) {
  const auto& loop_sequence = lut_bundle.LoopSequence();
  return TrimDuplex{lut_bundle.Sid5pMatcher(), lut_bundle.Sid3pMatcher(), lut_bundle.StartSequenceMatcher(),
                    lut_bundle.StopSequenceMatcher(), loop_sequence};
}

DemuxAndTrimDuplex::DemuxAndTrimDuplex(const LutBundleDuplex& lut_bundle) : _trim{CreateTrim(lut_bundle)} {}

void DemuxAndTrimDuplex::Demux(FixedReadRecord& record, DuplexMetrics& metrics) const {
  _trim.FindHairpin(record, metrics);
}

}  // namespace xoos::demux
