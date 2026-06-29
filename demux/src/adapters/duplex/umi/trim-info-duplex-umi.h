#pragma once

#include <optional>

#include "adapters/duplex/trim-info-duplex.h"

namespace xoos::demux {

/**
 * Extends TrimInfoDuplex with UMI (Unique Molecular Identifier) information.
 *
 * Example Structure:
 * CAACAATAGTTCAGACGTGTGCTCTTCCGATCT <umi-5p> <insert> <umi-3p> AGATCGGAAGAGCGTCGTGTAGG <hairpin>
 * CCTACACGACGCTCTTCCGATCT <rc-umi-3p> <rc-insert> <-rc-umi-5p> AGATCGGAAGAGCACACGTCTGAACTA
 */
struct TrimInfoDuplexUMI : TrimInfoDuplex {
  // UMI fields
  std::optional<uint> umi_5p;
  std::optional<uint> umi_3p;

  void Clear() {
    // Call parent class Clear() method
    TrimInfoDuplex::Clear();

    // Clear UMI fields
    umi_5p.reset();
    umi_3p.reset();
  }
};

}  // namespace xoos::demux
