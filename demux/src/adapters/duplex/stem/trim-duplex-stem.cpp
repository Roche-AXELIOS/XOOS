#include "trim-duplex-stem.h"

#include <xoos/log/logging.h>
#include <xoos/types/int.h>

#include "adapters/duplex/duplex-match.h"
#include "adapters/duplex/trim-duplex.h"
#include "io/read-record.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

/**
 * @brief Constructs a TrimDuplexStem adapter with matchers for duplex trimming and stem detection.
 *
 * @param sid_5p_matcher Matcher for the 5' sample identification sequence (relative to R1)
 * @param sid_3p_matcher Matcher for the 3' sample identification sequence (relative to R1)
 * @param start_matcher Matcher for the duplex start adapter sequence (relative to R1)
 * @param end_matcher Matcher for the duplex end adapter sequence (relative to R1)
 * @param stem_5p_matcher Matcher for the 5' stem sequence used to locate stem boundaries (relative to R1)
 * @param stem_3p_matcher Matcher for the 3' stem sequence used to locate stem boundaries (relative to R1)
 * @param loop_sequence The loop sequence separating 5' and 3' regions in the duplex
 */
TrimDuplexStem::TrimDuplexStem(SeqMatcher sid_5p_matcher, SeqMatcher sid_3p_matcher, SeqMatcher start_matcher,
                               SeqMatcher end_matcher, const SeqMatcher& stem_5p_matcher,
                               const SeqMatcher& stem_3p_matcher, std::string_view loop_sequence)
    : TrimDuplex(std::move(sid_5p_matcher), std::move(sid_3p_matcher), std::move(start_matcher), std::move(end_matcher),
                 loop_sequence),
      stem_5p_matcher(stem_5p_matcher),
      stem_3p_matcher(stem_3p_matcher),
      stem_5p_bitap(stem_5p_matcher.Pool().front().sequence, false),
      stem_3p_bitap(stem_3p_matcher.Pool().front().sequence, true) {}

/**
 * @brief Finds the 5' and 3' stem positions in the read.
 * Uses Bitap to find the position of the 5' and 3' stem sequences.
 * @param record a FixedReadRecord containing the sequence to search
 * @return a pair of integers representing the positions of the 5' and 3' stem to trim, and -1 if not found.
 */
std::pair<s32, s32> TrimDuplexStem::FindStem(FixedReadRecord& record) const {
  const auto seq = std::string_view(record.consensus_buffer, record.consensus_seq_len);
  // Find adapter sequence to anchor search using bitap
  // search using the 3' stem since it is more likely to exist
  s32 pos_3p = stem_3p_bitap.ReverseScan(seq);
  // if we found the 3' stem, search for the 5' stem upstream
  // right to left to compensate for artefacts (repeated stems)
  s32 pos_5p = stem_5p_bitap.ReverseScan(seq, -1, pos_3p);

  if ((pos_3p != -1 && pos_5p != -1) && pos_5p >= pos_3p) {
    // invalid stem positions
    return {-1, -1};
  }

  return {pos_5p, pos_3p};
}

}  // namespace xoos::demux
