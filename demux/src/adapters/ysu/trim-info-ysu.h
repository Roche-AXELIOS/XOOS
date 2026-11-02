#pragma once

#include <optional>

#include <xoos/types/int.h>

#include "sequence/loci-range.h"

namespace xoos::demux {
/**
 * Barcode Ids for matches if found.
 */
struct TrimInfoYsu {
  std::optional<u32> sid;

  std::optional<u32> sid_5p;
  std::optional<u32> sid_5p_edist;
  std::optional<u32> umi_5p;

  std::optional<u32> sid_3p;
  std::optional<u32> sid_3p_edist;
  std::optional<u32> umi_3p;

  LociRange insert;

  void Clear() {
    sid.reset();
    sid_5p.reset();
    sid_5p_edist.reset();
    umi_5p.reset();
    sid_3p.reset();
    sid_3p_edist.reset();
    umi_3p.reset();
    insert.Clear();
  }
};
}  // namespace xoos::demux
