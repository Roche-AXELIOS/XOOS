#pragma once

#include "io/read-record.h"
#include "lut-bundle-ysu-te.h"
#include "trim-info-ysu-te.h"
#include "trim3p-ysu-te.h"
#include "trim5p-ysu-te.h"

namespace xoos::demux {

/**
 * @class DemuxAndTrimYsuTe
 * @brief Trims 3' and 5' YSU-TE adapters, identifies SID and UMI, and determines insert start and end positions.
 *
 * This class is responsible for processing reads by trimming adapter sequences from both the 3' and 5' ends,
 * identifying Sample ID (SID) and Unique Molecular Identifier (UMI), and determining the start and end positions
 * of the insert sequence. The logic is specific to the YSU-TE version of adapter architecture.
 */
class DemuxAndTrimYsuTe {
 public:
  /**
   * Construct a DemuxAndTrimYsuTe object.
   *
   * @param enable_partial Enable partial reads (reads with partial or no 3' adapter)
   * @param lut_bundle The bundle of LUTs for all necessary LUTs
   */
  DemuxAndTrimYsuTe(bool enable_partial, const LutBundleYsuTe& lut_bundle);

  /**
   * Trims adapters, identifies SID and UMI, determines insert start and end, and determines which sample to
   * assign the read.
   *
   * @param record The read to be trimmed and assigned to a sample
   * @return A structure containing the determined 5' SID & UMI, 3' SID & UMI, insert start & end, and assigned sample
   */
  TrimInfoYsuTe operator()(const FixedReadRecord& record) const;

 private:
  /**
   * Given the results of trimming the 5' and 3' adapters, determine the sample to assign the record to, and produce
   * the final trim results.
   *
   * @param trim_5p 5' trim results
   * @param trim_3p 3' trim results
   * @return A structure containing the determined 5' SID & UMI, 3' SID & UMI, insert start & end, and assigned sample
   */
  TrimInfoYsuTe Demux(const Trim5pInfoYsuTe& trim_5p, const Trim3pInfoYsuTe& trim_3p) const;

  /**
   * Determine if a read should be assigned to a sample, and if so determine which sample. The following criteria is
   * used to determine if a read should be assigned to a sample:
   *
   *    full length reads: at least 1 SID identified and both UMI identified
   *    partial reads: at least 1 SID and 1 UMI identified
   *
   * @param trim_5p 5' trim results
   * @param trim_3p 3' trim results
   * @return Which sample, if any, the read should be assigned. If std::nullopt, then do not assign read to sample.
   */
  std::optional<uint> DetermineSampleId(const Trim5pInfoYsuTe& trim_5p, const Trim3pInfoYsuTe& trim_3p) const;

 private:
  bool _enable_partial;
  Trim5pYsuTe _trim_5p;
  Trim3pYsuTe _trim_3p;
};
}  // namespace xoos::demux
