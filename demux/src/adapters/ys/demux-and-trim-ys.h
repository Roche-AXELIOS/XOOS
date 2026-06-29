#pragma once
#include "io/read-record.h"
#include "lut-bundle-ys.h"
#include "trim-info-ys.h"
#include "trim3p-ys.h"
#include "trim5p-ys.h"

namespace xoos::demux {
/**
 * @class DemuxAndTrimYs
 * @brief Trims 3' and 5' YS adapters, identifies SID, spacer, and determines insert start and end positions.
 *
 * This class is responsible for processing reads by trimming adapter sequences from both the 3' and 5' ends,
 * identifying Sample ID (SID), and determining the start and end positions
 * of the insert sequence. The logic is specific to the YS version of adapter architecture.
 */
class DemuxAndTrimYs {
 public:
  /**
   * Construct a DemuxAndTrimYs object.
   *
   * @param enable_partial Enable partial reads (reads with partial or no 3' adapter)
   * @param lut_bundle The bundle of LUTs for all necessary LUTs
   */
  DemuxAndTrimYs(bool enable_partial, const LutBundleYs& lut_bundle);
  /**
   * Trims adapters, identifies SID, determines insert start and end, and determines which sample to
   * assign the read.
   *
   * @param record The read to be trimmed and assigned to a sample
   * @return A structure containing the determined 5' SID, 3' SID, insert start & end, and assigned sample
   */
  TrimInfoYs operator()(const FixedReadRecord& record) const;

 private:
  /**
   * Given the results of trimming the 5' and 3' adapters, determine the sample to assign the record to, and produce
   * the final trim results.
   *
   * @param trim_5p 5' trim results
   * @param trim_3p 3' trim results
   * @return A structure containing the determined 5' SID, 3' SID, insert start & end, and assigned sample
   */
  TrimInfoYs Demux(const Trim5pInfoYs& trim_5p, const Trim3pInfoYs& trim_3p) const;

  /**
   * Determine if a read should be assigned to a sample, and if so determine which sample. The following criteria is
   * used to determine if a read should be assigned to a sample:
   *
   *    full length reads: Both SID identified
   *
   *    partial reads: at least 1 SID (unlike read with UMI -> least 1 SID and 1 UMI identified)
   *
   * @param trim_5p 5' trim results
   * @param trim_3p 3' trim results
   * @return Which sample, if any, the read should be assigned. If std::nullopt, then do not assign read to sample.
   */
  std::optional<u32> DetermineSampleId(const Trim5pInfoYs& trim_5p, const Trim3pInfoYs& trim_3p) const;

  bool _enable_partial;
  Trim5pYs _trim_5p;
  Trim3pYs _trim_3p;
};
}  // namespace xoos::demux
