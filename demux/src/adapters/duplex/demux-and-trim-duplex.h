#pragma once

#include <xoos/types/int.h>

#include "io/read-record.h"
#include "lut-bundle-duplex.h"
#include "trim-duplex.h"

namespace xoos::demux {
struct DuplexMetrics;

/**
 * @class DemuxAndTrimDuplex
 * @brief Trims 3' and 5' v9.2 adapters, identifies SID and UMI, and determines insert start and end positions.
 *
 * This class is responsible for processing reads by trimming adapter sequences from both the 3' and 5' ends,
 * identifying Sample ID (SID) and Unique Molecular Identifier (UMI), and determining the start and end positions
 * of the insert sequence. The logic is specific to the v9.2 version of adapter architecture.
 */
class DemuxAndTrimDuplex {
 public:
  virtual ~DemuxAndTrimDuplex() = default;
  /**
   * Construct a DemuxAndTrimDuplex object.
   *
   * @param lut_bundle The bundle of LUTs for all necessary LUTs
   */
  explicit DemuxAndTrimDuplex(const LutBundleDuplex& lut_bundle);

  /**
   * Find the start adapter in the consensus of a pair-wise aligned read. If not successful, attempt to find
   * the start adapter.
   */
  s32 FindStartAdapterInConsensus(const FixedReadRecord& record) const {
    return _trim.FindStartAdapterInConsensus(record);
  }

  /**
   * Pure virtual function to find the UMI position in a read.
   *
   * This function was introduced as a workaround for compatibility with the aligner interface for SIMDNA.
   * TODO: Remove or refactor this method when the aligner no longer requires it especially if duplex UMI is refactored
   */
  virtual std::pair<s32, s32> FindUMIPos(FixedReadRecord& /*unused*/) const {
    // Is a dummy
    return {};
  }

 private:
  TrimDuplex _trim;
};
}  // namespace xoos::demux
