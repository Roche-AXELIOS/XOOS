#include "trim-duplex.h"

#include <xoos/enum/enum-util.h>
#include <xoos/types/int.h>

#include <map>
#include <utility>

#include "duplex-match.h"

namespace xoos::demux {

TrimDuplex::TrimDuplex(SeqMatcher sid_5p_matcher, SeqMatcher sid_3p_matcher, SeqMatcher sid_start_matcher,
                       SeqMatcher sid_end_matcher, std::string_view loop_sequence)
    : _sid_5p_matcher(std::move(sid_5p_matcher)),
      _sid_3p_matcher{std::move(sid_3p_matcher)},
      _sid_start_matcher{std::move(sid_start_matcher)},
      _sid_end_matcher{std::move(sid_end_matcher)},
      _loop_sequence{loop_sequence},
      _cascaded_luts(_sid_5p_matcher, _sid_3p_matcher, DuplexMatch::BarcodeType::kSID5p,
                     DuplexMatch::BarcodeType::kSID3p),
      _mask_sid5p(static_cast<int64_t>((1ul << (2 * _cascaded_luts.MaxLength(DuplexMatch::kSID5p))) - 1ul)),
      _mask_sid3p(static_cast<int64_t>((1ul << (2 * _cascaded_luts.MaxLength(DuplexMatch::kSID3p))) - 1ul)),
      _mask_plan_a_full{_mask_sid5p, _mask_sid3p, 0, 0, 0, 0, 0, 0},
      _kStart(_sid_start_matcher.Pool().front().sequence, false),
      // Abbreviated sequence
      _kShortStart(static_cast<std::string_view>(_sid_start_matcher.Pool().front().sequence).substr(7), false),
      _kHairpin(_loop_sequence, false)
// Abbreviated sequence
{
  // Initialize the bitap's using a std::map, then convert to a std::vector. Note that for performance reasons, most
  // member variables in the bitap are const, so we need to initialize them upon construction which makes the use of
  // a vector more difficult.
  std::map<u32, Bitap<4>> bitap_5p, bitap_3p;
  for (auto& sid : _sid_5p_matcher.Pool()) {
    // We need to search in reverse order for the 5' barcode.
    bitap_5p.try_emplace(sid.id, Bitap<4>(sid.sequence, true));
  }
  for (auto& sid : _sid_3p_matcher.Pool()) {
    bitap_3p.try_emplace(sid.id, Bitap<4>(sid.sequence, false));
  }
  // The map will sort the ids, so figure out the max id and create a vector.
  // We should expect that SID list is not empty
  const auto max_id = bitap_5p.rbegin()->first;

  _sid_5p_bitap.reserve(max_id + 1);
  _sid_3p_bitap.reserve(max_id + 1);
  const Bitap<4> dummy("", true);
  for (int i = 0; std::cmp_less_equal(i, max_id); ++i) {
    auto iter{bitap_5p.find(i)};
    if (iter != bitap_5p.end()) {
      _sid_5p_bitap.push_back(iter->second);
      // FInd corresponding entry in bitap_3p, should be present
      _sid_3p_bitap.push_back(bitap_3p.find(i)->second);
    } else {
      _sid_5p_bitap.push_back(dummy);
      _sid_3p_bitap.push_back(dummy);
    }
  }
  // the bitap objects are now accessible by id in O(1) time.
}

int TrimDuplex::FindStartAdapterInConsensus(FixedReadRecord& record) const {
  // TODO: I made these constants from the magic numbers in the original code. I'm not sure if they are correct.
  //               We may want to reimplement this entire function. It is difficult to read and understand.
  constexpr int kMaxStartAdapterLength{18};                         // max length of start adapter
  constexpr int kPartialStartAdapterLength{12};                     // search length of a partial start adapter
  constexpr int kSearchLength{63};                                  // length of the search window (64 - 1)
  constexpr int kIncr{kSearchLength + 1 - kMaxStartAdapterLength};  // increment after each search
  constexpr int kDoubleReadThresh{450};                             // Threshold to check if it is a double read

  const int max_pos{record.consensus_seq_len - 1};  // last position in the consensus sequence
  int pos = -1;

  // Determine if we have a double read. If the hairpin position is greater than our threshold it may be a double read.
  if (max_pos > kDoubleReadThresh) {
    // we expect the hairpin position to be around 200. If the hairpin is at a higher position, we might be dealing
    // with a concatenated read. In that case, start searching for the start adapter starting from the midadapter
    // rather than the beginning of the read.
    int end = max_pos;
    int start = end - kSearchLength;
    while (pos == -1 && end >= 0) {
      pos = _kStart.Find(record.consensus_buffer, std::max(0, start), end);
      start -= kIncr;
      end -= kIncr;
    }
    if (pos != -1) {
      return pos;
    }
  } else {
    int start = 0;
    int end = kSearchLength;
    // The usual recipe: find start adapter starting at begin
    while (pos == -1 && start <= max_pos) {
      pos = _kStart.Find(record.consensus_buffer, start, std::min(end, max_pos));
      start += kIncr;  // 18 is max length of start adapter
      end += kIncr;
    }
    if (pos != -1 && pos < max_pos) {
      return pos;
    }
  }

  // we did not find a start adapter. It's likely not there, but look for a truncated adapter at the beginning.
  return _kShortStart.Find(record.consensus_buffer, 0, std::min(kPartialStartAdapterLength, max_pos));
}
}  // namespace xoos::demux
