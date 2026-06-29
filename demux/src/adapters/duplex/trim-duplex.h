#pragma once

#include <xoos/types/int.h>

#include <array>
#include <vector>

#include "adapters/duplex/duplex-match.h"
#include "io/read-record.h"
#include "sequence/matcher/bitap.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

struct DuplexMetrics;

// TODO: hardcoded values should be removed, either taken from configuration or dynamically determined
// this should be based on the length of the sequence being searched for
constexpr s32 kBitapErrors = 4;

class TrimDuplex {
 public:
  TrimDuplex(SeqMatcher sid_5p_matcher, SeqMatcher sid_3p_matcher, SeqMatcher sid_start_matcher,
             SeqMatcher sid_end_matcher, std::string_view loop_sequence);

  ~TrimDuplex() = default;

  // This is the first step of Duplex trip: finding the (likely) center of the read, which should contain a hairpin.
  // Code implemented in hairpin.cpp. Results of this step are written into the record that is passed in.
  void FindHairpin(FixedReadRecord& record, DuplexMetrics& metrics) const;

  // Find the start adapter in the consensus read
  int FindStartAdapterInConsensus(FixedReadRecord& record) const;

 private:
  SeqMatcher _sid_5p_matcher;
  SeqMatcher _sid_3p_matcher;
  SeqMatcher _sid_start_matcher;
  SeqMatcher _sid_end_matcher;
  std::string_view _loop_sequence;
  const CascadedLUTs _cascaded_luts;
  const s64 _mask_sid5p;
  const s64 _mask_sid3p;
  // Various pairs of values used during the search for the hairpin (sids) and start/end adapters.
  // plan_a is the initial filter check
  const std::array<s64, 8> _mask_plan_a_full;
  const std::array<s64, 8> _types_plan_a_full = {
      DuplexMatch::BarcodeType::kSID5p,   DuplexMatch::BarcodeType::kSID3p,   DuplexMatch::BarcodeType::kUnknown,
      DuplexMatch::BarcodeType::kUnknown, DuplexMatch::BarcodeType::kUnknown, DuplexMatch::BarcodeType::kUnknown,
      DuplexMatch::BarcodeType::kUnknown, DuplexMatch::BarcodeType::kUnknown};

  // Determine position of midadapter "hairpin".
  void FilterResults(FixedReadRecord& record, const CascadedLUTs& cascaded_luts) const;

  // After application of the cascaded LUTs, we likely will have found SIDs - this function determines which
  // pairs of 5p and 3p give the lowest error (and should be the best match).
  int CalculateMinimumErrorHairpin(FixedReadRecord& record, const CascadedLUTs& cascaded_luts, int& best_match_index_5,
                                   int& best_match_index_3) const;

  // Encapsulation of Bitap alignment algorithm to find start and end markers.
  // Edit distance tolerances are set to at least an error rate > 20% and at most 30% (~25% is target).
  // TODO: consider making this based on probability of random matches rather than fixed error rates (requires testing)
  // 4/18 = 22% error rate as 18bp is the expect length of the start adapter
  const Bitap<4> _kStart;
  // 3/12 = 25% error rate but 12 is a length determined empirically to support truncated adapters
  // Only searched for if the full-length start adapter is not found.
  const Bitap<3> _kShortStart;
  // 2/7 = 28% error rate as 7bp is the expect length of the loop sequence
  // Note that the loop is only 7bp long we later use the detection of SIDs to improve specificity.
  const Bitap<2> _kHairpin;

  // We also can use Bitap to find the SID adapters. As the bitap algorithm finds the "end position" of the
  // match, we will search in the reverse direction for the 5p adapter and in the forward direction for the 3p adapter.
  std::vector<Bitap<4>> _sid_5p_bitap;
  std::vector<Bitap<4>> _sid_3p_bitap;

  // Helper struct to pass around some results of the duplex trim.
  struct TrimResults {
    u32 length{0};
  };

  // First step of duplex trim: look for the kmer "ATGTATG" using a string search and look for SIDs next to it.
  void FindHairpinByStringSearch(FixedReadRecord& record) const;

  // Second step of the duplex trim: after finding global symmetry in the read (very fast), look for hairpin
  // and SID adapters after that.
  void FindHairpinByGlobalSymmetry(FixedReadRecord& record) const;

  // We're doing an extensive search for the hairpin in the read  and look for the SID adapters if
  // we found the middle. T
  void FindHairpinByLocalSymmetry(FixedReadRecord& record, TrimResults& results) const;
};
}  // namespace xoos::demux
