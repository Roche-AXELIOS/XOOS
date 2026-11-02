#pragma once

#include <algorithm>
#include <span>
#include <string_view>

#include <xoos/types/int.h>

namespace xoos::demux {

/**
 * @brief Indel Smearing: pseudo-affine gap topology optimization for Levenshtein alignments.
 *
 * Edlib produces unit-cost (Levenshtein) pairwise alignments where insertions, deletions, and
 * mismatches each cost 1. This means scattered single-base indels score identically to grouped
 * multi-base gaps — the alignment has no preference for gap grouping.
 *
 * Indel smearing addresses this by post-processing the alignment array: indels are slid along
 * the consensus sequence in a chosen direction without changing the edit distance or consensus
 * bases. When two indels collide (slide into each other), they merge into a single longer gap,
 * reducing the total number of gap openings. This produces a gap topology closer to what an
 * affine-gap-penalty aligner would yield.
 *
 * Directionality:
 *   - SmearIndelsLeft:  iterates right-to-left, pushing indels toward index 0.
 *   - SmearIndelsRight: iterates left-to-right, pushing indels toward the end.
 *   - Group smearing (SmearIndelsLeft then SmearIndelsRight, or vice versa): first collapses
 *     scattered indels together, then repositions the group to the desired strand-relative anchor.
 *
 * For strand-aware processing (Alignment::SmearAlignment):
 *   - Prefix mode produces a reversed alignment array. The pass order swaps relative to infix
 *     mode so that grouped indels land on the same anchor side relative to the read structure.
 *   - When strand detection is enabled, forward and reverse strands get opposite final positioning
 *     so that both strands' indels land consistently relative to the read.
 *
 * Candidate selection (Alignment::SelectBestPrefixAlign):
 *   When multiple equal-edit-distance tracebacks exist, the candidate with the fewest total
 *   indel operations is selected. Since all candidates share the same edit distance, fewer
 *   indels means more mismatches — a cleaner topology for downstream quality scoring.
 *
 * Future:
 *   Indel smearing can be used to detect bases that have ambiguous alignment positions. If we
 *   smear left then right we can flag all bases that were within the smear as ambiguous even if
 *   strand detection is not active, which we can assign a different base quality.
 */

/**
 * @brief Runtime-configurable quality scores (phred+33 chars) for each base category.
 *
 * Defaults reproduce the historical hardcoded scores. The simplex score is applied on both the
 * simplex-adapter output path and the simplex-tail region of duplex reads; the other three apply
 * only to the duplex consensus. Overridable via --quality-score-overrides.
 */
struct BaseQualityScores {
  // indexes for each base replacement
  static constexpr auto kSimplexIndex = 0;
  static constexpr auto kConcordantIndex = 1;
  static constexpr auto kDiscordantIndex = 2;
  static constexpr auto kAmbiguousIndex = 3;
  static constexpr auto kQualityBaseValues = 4;

  // quality base defaults
  static constexpr char kDiscordBaseQual{'&'};
  static constexpr char kSimplexBaseQual{'7'};
  static constexpr char kConcordBaseQual{'H'};
  static constexpr char kAmbiguousBaseQual{'9'};

  char simplex{kSimplexBaseQual};
  char concordant{kConcordBaseQual};
  char discordant{kDiscordBaseQual};
  char ambiguous{kAmbiguousBaseQual};
};

// Edit operations. Fixed underlying type (u8) so comparisons and assignments against the u8
// alignment buffers do not trip mixed enum/integer comparison compiler issues (MISRA M23_061).
enum : u8 {
  kEdopMatch = 0,
  kEdopInsert = 1,
  kEdopDelete = 2,
  kEdopMismatch = 3,
  kEdopMatchMethyl = 4
};

// YC tag components
constexpr char kDefaultJunction{'+'};
constexpr char kSwappedJunction{'-'};

// IUPAC-like base encoding table for duplex consensus. This table is used to encode two bases into one IUPAC-like base
// and is used to create a compact, lossless compression of regions R1 and R2 used in the alignment process. The
// IUPAC-like encoding describes the differences between each of the two reads and the consensus sequence and is
// therefore a 5x5 table because we also need to encode indels.
//                                    A    C    G    T    -
constexpr char kIUPACTable[5][5] = {{'*', 'M', 'R', 'W', 'I'},   // A
                                    {'B', '*', 'S', 'Y', 'P'},   // C
                                    {'D', 'V', '*', 'K', 'J'},   // G
                                    {'H', 'E', 'F', '*', 'X'},   // T
                                    {'L', 'Q', 'O', 'Z', '*'}};  // -

/// @brief Index in kIUPACTable for indels, which is the last row/column
constexpr u8 kIndelIndex = 4;

/**
 * @brief Smears indels leftward in the alignment array (toward index 0).
 *
 * Iterates right-to-left, sliding each indel as far left as possible without introducing a
 * mismatch in the consensus. When indels collide they merge, reducing gap openings.
 * Does not modify the consensus sequence — only repositions operations in the alignment array.
 *
 * @param[in]     consensus The consensus sequence (duplex region only).
 * @param[in,out] alignment The alignment operation array to modify. Must be at least as long as consensus.
 */
void SmearIndelsLeft(std::string_view consensus, std::span<u8> alignment);

/**
 * @brief Smears indels rightward in the alignment array (toward the end).
 *
 * Iterates left-to-right, sliding each indel as far right as possible without introducing a
 * mismatch in the consensus. When indels collide they merge, reducing gap openings.
 * Does not modify the consensus sequence — only repositions operations in the alignment array.
 *
 * @param[in]     consensus The consensus sequence (duplex region only).
 * @param[in,out] alignment The alignment operation array to modify. Must be at least as long as consensus.
 */
void SmearIndelsRight(std::string_view consensus, std::span<u8> alignment);

/**
 * @brief Smears indels leftward and marks every position crossed (any repeat, incl. tandem).
 *
 * Use when indels were right-anchored and you want to find how far left they could slide.
 * The alignment span is mutated as a side effect (same as SmearIndelsLeft). The swap gate
 * (consensus[pos] == consensus[indel_pos]) makes the slide consensus-preserving; any repeat that
 * satisfies it is marked.
 *
 * @param[in]     consensus The consensus sequence (duplex region only).
 * @param[in,out] alignment The alignment operation array; will be mutated by the sweep.
 * @param[out]    mask      Output mask; 1 = ambiguous, 0 = not. Must be at least alignment.size().
 */
inline void SmearMarkAllLeft(const std::string_view consensus,
                             const std::span<u8> alignment,
                             const std::span<u8> mask) {
  const auto len = alignment.size();
  std::fill_n(mask.data(), len, u8{0});
  size_t indel_length = 0;
  for (auto consensus_pos = static_cast<s32>(len) - 1; consensus_pos >= 0; --consensus_pos) {
    bool is_indel = false;
    switch (alignment[static_cast<u32>(consensus_pos)]) {
      case kEdopInsert:
      case kEdopDelete:
        ++indel_length;
        is_indel = true;
        break;
      default:
        break;
    }
    if (!is_indel && indel_length > 0) {
      const auto indel_pos = static_cast<size_t>(consensus_pos) + indel_length;
      if (consensus[static_cast<size_t>(consensus_pos)] == consensus[indel_pos]) {
        mask[static_cast<size_t>(consensus_pos)] = 1;
        mask[indel_pos] = 1;
        std::swap(alignment[indel_pos], alignment[static_cast<size_t>(consensus_pos)]);
      } else {
        indel_length = 0;
      }
    }
  }
}

/**
 * @brief Smears indels rightward and marks every position crossed (any repeat, incl. tandem).
 *
 * Use when indels were left-anchored and you want to find how far right they could slide.
 * The alignment span is mutated as a side effect (same as SmearIndelsRight). The swap gate
 * (consensus[pos] == consensus[indel_pos]) makes the slide consensus-preserving; any repeat that
 * satisfies it is marked.
 *
 * @param[in]     consensus The consensus sequence (duplex region only).
 * @param[in,out] alignment The alignment operation array; will be mutated by the sweep.
 * @param[out]    mask      Output mask; 1 = ambiguous, 0 = not. Must be at least alignment.size().
 */
inline void SmearMarkAllRight(const std::string_view consensus,
                              const std::span<u8> alignment,
                              const std::span<u8> mask) {
  const auto len = alignment.size();
  std::fill_n(mask.data(), len, u8{0});
  u32 indel_length = 0;
  for (u32 consensus_pos = 0; consensus_pos < len; ++consensus_pos) {
    bool is_indel = false;
    switch (alignment[consensus_pos]) {
      case kEdopInsert:
      case kEdopDelete:
        ++indel_length;
        is_indel = true;
        break;
      default:
        break;
    }
    if (!is_indel && indel_length > 0) {
      const u32 indel_pos = consensus_pos - indel_length;
      if (consensus[consensus_pos] == consensus[indel_pos]) {
        mask[consensus_pos] = 1;
        mask[indel_pos] = 1;
        std::swap(alignment[indel_pos], alignment[consensus_pos]);
      } else {
        indel_length = 0;
      }
    }
  }
}

void MethylConvert(
    char* consensus, u8* alignment, s32 alignment_length, char* xm_tag_start, const u8* r2, bool swapped);

};  // namespace xoos::demux
