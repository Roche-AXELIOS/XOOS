#pragma once

#include <optional>

#include "sequence/loci-range.h"

namespace xoos::demux {
/**
 * Barcode Ids for matches if found.
 */
struct TrimInfoYs {
  std::optional<u32> sid;

  std::optional<u32> sid_5p;
  std::optional<u32> sid_5p_edist;

  std::optional<u32> sid_3p;
  std::optional<u32> sid_3p_edist;

  LociRange insert;

  void Clear() {
    sid.reset();
    sid_5p.reset();
    sid_5p_edist.reset();
    sid_3p.reset();
    sid_3p_edist.reset();
    insert.Clear();
  }
};
}  // namespace xoos::demux
