#pragma once

#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

namespace fs = std::filesystem;

class LutBundleDuplex {
 public:
  LutBundleDuplex(SeqMatcher sid_5p, SeqMatcher start_sequence, SeqMatcher sid_3p, SeqMatcher stop_sequence,
                  SeqMatcher loop_sequence);

  const SeqMatcher& Sid5pMatcher() const;
  const BarcodePool& Sid5pPool() const;

  const SeqMatcher& StartSequenceMatcher() const;

  const SeqMatcher& Sid3pMatcher() const;

  const SeqMatcher& StopSequenceMatcher() const;

  std::string_view LoopSequence() const;

 private:
  SeqMatcher _sid_5p_matcher;
  SeqMatcher _start_sequence_matcher;
  SeqMatcher _sid_3p_matcher;
  SeqMatcher _stop_sequence_matcher;
  SeqMatcher _loop_sequence_matcher;
};
}  // namespace xoos::demux
