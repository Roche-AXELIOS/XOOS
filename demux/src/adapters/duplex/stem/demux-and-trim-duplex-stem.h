#pragma once

#include "adapters/duplex/demux-and-trim-duplex.h"
#include "lut-bundle-duplex-stem.h"
#include "trim-duplex-stem.h"

namespace xoos::demux {
struct DuplexMetrics;

class DemuxAndTrimDuplexStem : public DemuxAndTrimDuplex {
 public:
  explicit DemuxAndTrimDuplexStem(const LutBundleDuplexStem& lut_bundle);

  std::pair<s32, s32> FindUMIPos(FixedReadRecord& read) const override;

 private:
  TrimDuplexStem _trim_stem;
};

}  // namespace xoos::demux
