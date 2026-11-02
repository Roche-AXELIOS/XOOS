#pragma once

#include <xoos/types/int.h>

#include "io/read-record.h"
#include "lut-bundle-ysu.h"
#include "trim-info-ysu.h"
#include "trim3p-ysu.h"
#include "trim5p-ysu.h"

namespace xoos::demux {

/**
 * @class DemuxAndTrimYsu
 * @brief Trims 3' and 5' YSU adapters, identifies SID and UMI, and determines insert start and end positions.
 *
 * This class is responsible for processing reads by trimming adapter sequences from both the 3' and 5' ends,
 * identifying Sample ID (SID) and Unique Molecular Identifier (UMI), and determining the start and end positions
 * of the insert sequence. The logic is specific to the YSU version of adapter architecture.
 */
class DemuxAndTrimYsu {
 public:
  /**
   * Construct a DemuxAndTrimYsu object.
   *
   * @param lut_bundle The bundle of LUTs for all necessary LUTs
   */
  explicit DemuxAndTrimYsu(const LutBundleYsu& lut_bundle);

  /**
   * Trims adapters, identifies SID and UMI, determines insert start and end, and determines which sample to
   * assign the read.
   *
   * @param record The read to be trimmed and assigned to a sample
   * @return A structure containing the determined 5' SID & UMI, 3' SID & UMI, insert start & end, and assigned sample
   */
  TrimInfoYsu operator()(const FixedReadRecord& record) const;

 private:
  /**
   * Given the results of trimming the 5' and 3' adapters, determine the sample to assign the record to, and produce
   * the final trim results. When SIDs are discordant, only the winning side is trimmed.
   *
   * @param trim_5p 5' trim results
   * @param trim_3p 3' trim results
   * @param seq_length Raw sequence length (used as untrimmed boundary for losing side)
   * @return A structure containing the determined 5' SID & UMI, 3' SID & UMI, insert start & end, and assigned sample
   */
  static TrimInfoYsu Demux(const Trim5pInfoYsu& trim_5p, const Trim3pInfoYsu& trim_3p, u32 seq_length);

  /**
   * Determine if a read should be assigned to a sample, and if so determine which sample. Assignment is SID-only: a
   * read is assigned whenever at least one SID is identified. The UMI requirement is enforced downstream by the
   * post-assignment umi_failed gate.
   *
   * @param trim_5p 5' trim results
   * @param trim_3p 3' trim results
   * @return Which sample, if any, the read should be assigned. If std::nullopt, then do not assign read to sample.
   */
  static std::optional<u32> DetermineSampleId(const Trim5pInfoYsu& trim_5p, const Trim3pInfoYsu& trim_3p);

  Trim5pYsu _trim_5p;
  Trim3pYsu _trim_3p;
};
}  // namespace xoos::demux
