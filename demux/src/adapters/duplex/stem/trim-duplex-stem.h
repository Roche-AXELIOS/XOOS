#pragma once

#include "adapters/duplex/trim-duplex.h"
#include "io/read-record.h"

namespace xoos::demux {

class TrimDuplexStem : public TrimDuplex {
 public:
  TrimDuplexStem(SeqMatcher sid_5p_matcher, SeqMatcher sid_3p_matcher, SeqMatcher start_matcher, SeqMatcher end_matcher,
                 const SeqMatcher& stem_5p_matcher, const SeqMatcher& stem_3p_matcher, std::string_view loop_sequence);

  std::pair<s32, s32> FindStem(FixedReadRecord& record) const;

 protected:
  const SeqMatcher& stem_5p_matcher;
  const SeqMatcher& stem_3p_matcher;

  // functionally the start adapter (minus runway)
  const Bitap<kBitapErrors> stem_5p_bitap;

  // novel sequence before loop (with parameter reverse = true so we get the start position of the match)
  const Bitap<kBitapErrors> stem_3p_bitap;
};

}  // namespace xoos::demux
