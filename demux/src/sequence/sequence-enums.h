#pragma once

namespace xoos::demux {

enum class ReadEnd {
  k5p,
  k3p,
};

enum class SeqTransform {
  kReverse,
  kReverseComplement,
  kNone,
};
}  // namespace xoos::demux
