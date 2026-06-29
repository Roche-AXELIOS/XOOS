#pragma once

#include <array>
#include <optional>

#include "adapters/duplex/duplex-match.h"
#include "sequence/loci-range.h"

namespace xoos::demux {
/**
 * Barcode Ids for matches if found. The recipe for duplex info is as follows:
 * {1 if mt_start_adaptor_found else 0}:{1 if mt_midadapter_found else 0}:{1 if mt_end_adapter_found else 0}:
 * {(mt_midadapter_start - mt_start_adapter_end) if (mt_midadapter_found and mt_start_adapter_found) else
 * (mt_shorted_read1_end - mt_start_adapter_end)
 * if mt_shorted_read_found else -1}:{(mt_midadapter_end - mt_start_adapter_end)
 * if (mt_midadapter_found and mt_start_adapter_found) else (mt_shorted_read2_start - mt_start_adapter_end)
 * if mt_shorted_read_found else -1}",
 */
struct TrimInfoDuplex  // NOLINT - unfiltered_matches does not need to be initialized
{
  // Status of duplex trimming is determined in hairpin.cpp.
  // kNone: Not enough markers found to do anything useful (notably no start and 5p SID marker found)
  // kZeroPlus: Found "half" of the markers (start and 5p SID marker/3p SID marker), but not enough 3p data
  enum class DuplexStatus { kNone, kZeroPlus, kMidAdapterFound };
  DuplexStatus duplex_status = DuplexStatus::kNone;

  // Used for Duplex HD trimming
  std::array<DuplexResult, DuplexMatch::kUnknown>
      matches;  // For each duplex match type - "unknown" is # of match types
  uint symmetry_5p_pos = 0;
  uint symmetry_3p_pos = 0;
  uint symmetry_middle = 0;

  std::optional<LociRange> midadapter_range;

  uint32_t unfiltered_matches_count = 0;
  // Should be far more than enough for any reasonable case
  static const int kMaxUnfilteredMatches = 512;
  std::array<DuplexResult, kMaxUnfilteredMatches> unfiltered_matches;
  static const int kMaxFilteredMatches = 128;
  std::array<std::array<DuplexResult, kMaxFilteredMatches>, DuplexMatch::kUnknown> filtered_matches;
  std::array<uint32_t, DuplexMatch::kUnknown> filtered_matches_count;

  void Clear() {
    midadapter_range.reset();
    for (auto& match : matches) {
      match = DuplexResult{};
    }
    duplex_status = DuplexStatus::kNone;
    symmetry_5p_pos = symmetry_3p_pos = symmetry_middle = 0;
    unfiltered_matches_count = 0;
  }
};
}  // namespace xoos::demux
