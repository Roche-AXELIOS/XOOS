#pragma once

#include "adapters/duplex/trim-duplex.h"
#include "io/read-record.h"
#include "sequence/fixed-sequence.h"

namespace xoos::demux {

class TrimDuplexStem : public TrimDuplex {
 public:
  TrimDuplexStem(FixedSequence runway, FixedSequence start, const FixedSequence& stem);

  std::pair<s32, s32> FindStem(FixedReadRecord& record) const;

 protected:
  const FixedSequence& Stem() const;
  const Bitap<kBitapErrors>& StartBitap() const;
  const Bitap<kBitapErrors>& StemBitap() const;

 private:
  FixedSequence _stem;

  // start/anchor sequence immediately before the insert
  const Bitap<kBitapErrors> _start_bitap;

  // novel sequence right before hairpin
  const Bitap<kBitapErrors> _stem_bitap;
};

}  // namespace xoos::demux
