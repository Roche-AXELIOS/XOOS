#pragma once

#include "adapters/duplex/stem/demux-and-trim-duplex-stem.h"
#include "lut-bundle-duplex-umi.h"
#include "trim-duplex-umi.h"

namespace xoos::demux {
struct DuplexMetrics;

class DemuxAndTrimDuplexUmi : public DemuxAndTrimDuplexStem {
 public:
  explicit DemuxAndTrimDuplexUmi(const LutBundleDuplexUmi& lut_bundle);

  std::pair<s32, s32> FindUMIPos(FixedReadRecord& read) const override;

 private:
  TrimDuplexUmi _trim_umi;
};

}  // namespace xoos::demux
