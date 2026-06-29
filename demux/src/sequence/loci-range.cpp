#include "loci-range.h"

#include "utility/string-util.h"

namespace xoos::demux {
LociRange::LociRange() : spos{0}, epos{0} {}

LociRange::LociRange(uint spos, uint epos) : spos(spos), epos(epos) {}

uint LociRange::Length() const { return epos - spos; }

std::string LociRange::Sequence(const std::string& sequence) const { return std::string(Substr(sequence, spos, epos)); }
}  // namespace xoos::demux
