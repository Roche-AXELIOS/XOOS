#pragma once

#include "lut-bundle/lut-bundle.h"
#include "sequence/fixed-sequence.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * LUT bundle for simplex-10x adapter (5'-only).
 *
 * Read structure: [fixed_1 (runway)][sid_2 (SID)][fixed_3 (stem)][insert...]
 * No 3' adapter components. The stem (fixed_3 / p_read_1t) is used for
 * scoring only — the trim position is set after the SID, not after the stem.
 */
class LutBundleSimplex10x {
 public:
  LutBundleSimplex10x(FixedSequence fixed_1, SeqMatcher sid_2, FixedSequence fixed_3);
  LutBundleSimplex10x(LutBundleSimplex10x&&) noexcept = default;
  LutBundleSimplex10x(const LutBundleSimplex10x&) = delete;
  LutBundleSimplex10x& operator=(LutBundleSimplex10x&&) = delete;
  LutBundleSimplex10x& operator=(const LutBundleSimplex10x&) = delete;

  // 1: anchor / runway (e.g. CAACAA)
  FixedSequence fixed_1;
  // 2: SID barcode (12 bp)
  SeqMatcher sid_2_matcher;
  // 3: stem / p_read_1t (e.g. CTACACGACACTCTTCCGATCT, 22 bp) — scoring only, not trimmed.
  FixedSequence fixed_3;
};

template <>
LutBundleSimplex10x CreateLutBundle(const AdapterDesignBundle& designs);  // NOSONAR(M23_184)

}  // namespace xoos::demux
