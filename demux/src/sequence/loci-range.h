#pragma once

#include <xoos/types/int.h>

#include <string>

namespace xoos::demux {
/**
 * Represent a range of loci, the range is closed
 * at the start position, and open at the end position.
 * In other words: [spos, epos)
 */
struct LociRange {
  u32 spos;
  u32 epos;

  LociRange();
  LociRange(u32 spos, u32 epos);

  u32 Length() const;

  std::string Sequence(const std::string& sequence) const;

  void Clear() { spos = epos = 0; }
};

}  // namespace xoos::demux
