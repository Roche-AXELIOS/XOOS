#include "trim-duplex-stem.h"

#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "adapters/duplex/trim-duplex.h"
#include "io/read-record.h"
#include "sequence/matcher/match-position-states.h"

namespace xoos::demux {

/**
 * @brief Constructs a TrimDuplexStem adapter with fixed sequences for duplex trimming and stem detection.
 *
 * @param runway The runway sequence (relative to R1)
 * @param start The start/anchor sequence (relative to R1)
 * @param stem The sequence before the hairpin loop (relative to R1)
 */
// Positional constructor; adjacent FixedSequence params encode fixed read-structure positions
// documented in the class header and are constructed positionally by the assembler.
TrimDuplexStem::TrimDuplexStem(FixedSequence runway,  // NOSONAR(S5419)
                               FixedSequence start,
                               const FixedSequence& stem)
    : TrimDuplex(std::move(runway), std::move(start)),
      _stem(stem),
      _start_bitap(Start(), SearchDirection::kForward),
      _stem_bitap(_stem.Seq(), SearchDirection::kBackward) {
}

const FixedSequence& TrimDuplexStem::Stem() const {
  return _stem;
}

const Bitap<kBitapErrors>& TrimDuplexStem::StartBitap() const {
  return _start_bitap;
}

const Bitap<kBitapErrors>& TrimDuplexStem::StemBitap() const {
  return _stem_bitap;
}

/**
 * @brief Finds the start and stem anchor positions in the read.
 * Uses Bitap to locate the stem (via StemBitap) and the start/anchor sequence (via StartBitap).
 * @param record a FixedReadRecord containing the sequence to search
 * @return a pair of integers {start, end} representing the anchor positions to trim,
 *         with kNoMatchPosition for either value if not found or if the positions are invalid.
 */
std::pair<s32, s32> TrimDuplexStem::FindStem(FixedReadRecord& record) const {
  const auto seq = std::string_view(record.consensus_buffer, record.consensus_seq_len);
  // Find adapter sequence to anchor search using bitap
  s32 end = _stem_bitap.ReverseScan(seq, kNoMatchPosition, kNoMatchPosition);
  // right to left to compensate for artefacts (repeated starts)
  s32 start = _start_bitap.ReverseScan(seq, kNoMatchPosition, end);

  if ((start != kNoMatchPosition && end != kNoMatchPosition) && start >= end) {
    // invalid stem positions
    return {kNoMatchPosition, kNoMatchPosition};
  }

  return {start, end};
}

}  // namespace xoos::demux
