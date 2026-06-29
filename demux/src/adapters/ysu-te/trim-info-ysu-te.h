#pragma once

#include <optional>

#include "sequence/loci-range.h"

namespace xoos::demux {
/**
 * Barcode Ids for matches if found.
 */
struct TrimInfoYsuTe {
  std::optional<uint> sid;

  std::optional<uint> sid_5p;
  std::optional<uint> sid_5p_edist;
  std::optional<uint> umi_5p;

  std::optional<uint> sid_3p;
  std::optional<uint> sid_3p_edist;
  std::optional<uint> umi_3p;

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
