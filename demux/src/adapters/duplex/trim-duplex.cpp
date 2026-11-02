#include "trim-duplex.h"

#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include "sequence/matcher/match-position-states.h"

namespace xoos::demux {

constexpr size_t kShortStartAdapterOffset = 7;

TrimDuplex::TrimDuplex(FixedSequence runway, FixedSequence start)
    : _runway{std::move(runway)},
      _start{std::move(start)},
      _kStart(_start.Seq(), SearchDirection::kForward),
      // Abbreviated sequence
      _kShortStart(_start.Seq().substr(kShortStartAdapterOffset), SearchDirection::kForward) {
}

s32 TrimDuplex::FindStartAdapterInConsensus(const FixedReadRecord& record) const {
  // TODO: I made these constants from the magic numbers in the original code. I'm not sure if they are correct.
  //               We may want to reimplement this entire function. It is difficult to read and understand.

  // max length of start adapter
  const u8 max_start_adapter_length = static_cast<u8>(_start.Size()) + 2;
  // search length of a partial start adapter
  constexpr u8 kPartialStartAdapterSearchLength = 12;
  // length of the search window (64 - 1)
  constexpr u8 kConsensusSearchWindowLength = 63;
  // increment after each search
  const u16 search_increment = kConsensusSearchWindowLength + 1 - max_start_adapter_length;
  // Threshold to check if it is a double read
  constexpr u16 kDoubleReadHeuristicThreshold = 450;
  // last position in the consensus sequence
  const s32 max_pos{record.consensus_seq_len - 1};
  s32 pos = kNoMatchPosition;

  // Determine if we have a double read. If the hairpin position is greater than our threshold it may be a double read.
  if (max_pos > kDoubleReadHeuristicThreshold) {
    // we expect the hairpin position to be around 200. If the hairpin is at a higher position, we might be dealing
    // with a concatenated read. In that case, start searching for the start adapter starting from the midadapter
    // rather than the beginning of the read.
    s32 end = max_pos;
    s32 start = end - kConsensusSearchWindowLength;
    while (pos == kNoMatchPosition && end >= 0) {
      pos = _kStart.Find(record.consensus_buffer, std::max(0, start), end);
      start -= search_increment;
      end -= search_increment;
    }
    if (pos != kNoMatchPosition) {
      return pos;
    }
  } else {
    s32 start = 0;
    s32 end = kConsensusSearchWindowLength;
    // The usual recipe: find start adapter starting at begin
    while (pos == kNoMatchPosition && start <= max_pos) {
      pos = _kStart.Find(record.consensus_buffer, start, std::min(end, max_pos));
      start += search_increment;
      end += search_increment;
    }
    if (pos != kNoMatchPosition && pos < max_pos) {
      return pos;
    }
  }

  // we did not find a start adapter. It's likely not there, but look for a truncated adapter at the beginning.
  return _kShortStart.Find(record.consensus_buffer, 0, std::min<s32>(kPartialStartAdapterSearchLength, max_pos));
}

std::string_view TrimDuplex::RunwaySequence() const {
  return _runway.Seq();
}

std::string_view TrimDuplex::Start() const {
  return _start.Seq();
}

}  // namespace xoos::demux
