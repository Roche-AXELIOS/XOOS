#pragma once

#include <array>
#include <limits>
#include <optional>
#include <string_view>
#include <vector>

#include <xoos/types/int.h>

#include "adapters/duplex/duplex-match.h"
#include "adapters/duplex/lut-bundle-duplex.h"
#include "io/read-record.h"
#include "sequence/matcher/bitap.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

struct DuplexMetrics;

class HairpinFinder {
 public:
  HairpinFinder(SeqMatcher sid_5p_matcher, SeqMatcher sid_3p_matcher, std::string_view loop_sequence);

  /**
   * @brief Convenience constructor that pulls all matchers from a duplex LUT bundle.
   *
   * Also accepts LutBundleDuplexStem and LutBundleDuplexUmi via the LutBundleDuplex base reference,
   * which is why the stem/UMI adapters can share this single overload.
   */
  explicit HairpinFinder(const LutBundleDuplex& lut_bundle);

  ~HairpinFinder() = default;

  // This is the first step of Duplex trip: finding the (likely) center of the read, which should contain a hairpin.
  // Code implemented in hairpin.cpp. Results of this step are written into the record that is passed in.
  void FindHairpin(FixedReadRecord& record, DuplexMetrics& metrics) const;

 private:
  SeqMatcher _sid_5p_matcher;
  SeqMatcher _sid_3p_matcher;
  std::string_view _loop_sequence;
  const CascadedLUTs _cascaded_luts;
  const s64 _mask_sid5p;
  const s64 _mask_sid3p;
  // Various pairs of values used during the search for the hairpin (sids) and start/end adapters.
  // plan_a is the initial filter check
  const std::array<s64, 8> _mask_plan_a_full;
  const std::array<s64, 8> _types_plan_a_full = {DuplexMatch::BarcodeType::kSID5p,
                                                 DuplexMatch::BarcodeType::kSID3p,
                                                 DuplexMatch::BarcodeType::kSID5p,
                                                 DuplexMatch::BarcodeType::kSID3p,
                                                 DuplexMatch::BarcodeType::kUnknown,
                                                 DuplexMatch::BarcodeType::kUnknown,
                                                 DuplexMatch::BarcodeType::kUnknown,
                                                 DuplexMatch::BarcodeType::kUnknown};

  // Determine position of midadapter "hairpin".
  void FilterResults(FixedReadRecord& record, const CascadedLUTs& cascaded_luts) const;

  // After application of the cascaded LUTs, we likely will have found SIDs - this function determines which
  // pairs of 5p and 3p give the lowest error (and should be the best match).
  s32 CalculateMinimumErrorHairpin(FixedReadRecord& record,
                                   const CascadedLUTs& cascaded_luts,
                                   s32& best_match_index_5,
                                   s32& best_match_index_3) const;

  /**
   * @brief Aggregated result of the 5p/3p SID pair scoring pass.
   */
  struct SidPairScore {
    s32 min_error{std::numeric_limits<s32>::max()};
    s32 error5p{std::numeric_limits<s32>::max()};
    u32 min_5p_index{0};
    s32 best_5p{kNoMatchPosition};
    s32 best_3p{kNoMatchPosition};
  };

  /**
   * @brief Scores all 5p/3p SID candidate pairs by proximity to the hairpin loop and edit distance.
   *
   * Iterates over filtered 5p matches, pairs each with compatible 3p matches (same barcode_id),
   * and tracks the combination with the lowest total error.
   *
   * @param[in]  trim_info     Filtered match data from hairpin detection.
   * @param[in]  cascaded_luts LUT bundle for barcode length lookups.
   * @param[in]  loop_end_pos  Inclusive end position of the hairpin loop in the read.
   * @param[in]  loop_length   Nominal loop sequence length.
   * @return Scoring result containing minimum error, best 5p/3p indices, and fallback state.
   */
  static SidPairScore ScoreSidPairs(const TrimInfoDuplex& trim_info,
                                    const CascadedLUTs& cascaded_luts,
                                    s32 loop_end_pos,
                                    s32 loop_length);

  /**
   * @brief Attempts bitap-based 3p SID recovery when the LUT found no 3p match.
   *
   * Only applies when no 3p candidates exist and the best 5p error is within the recoverable
   * threshold. Uses the per-SID forward+backward bitap pair keyed on the 5p's barcode_id.
   *
   * @param[in,out] record            Read record to populate with the recovered 3p match.
   * @param[in]     score             Result from ScoreSidPairs (provides min_5p_index and error5p).
   * @param[out]    best_match_index_3 Set to 0 on success.
   * @return The error value on success/failure, or std::nullopt if the fallback doesn't apply.
   */
  std::optional<s32> TryBitapFallback3p(FixedReadRecord& record,
                                        const SidPairScore& score,
                                        s32& best_match_index_3) const;

  // Bitap matcher for the loop sequence (edit distance ≤ 2). The loop is short, so we rely on
  // subsequent SID detection to improve specificity.
  const Bitap<2> _loop_fw;
  // Reverse direction counterpart of _loop_fw, used by FindStartEnd to recover the loop start position.
  const Bitap<2> _loop_bw;

  // Precomputed adapter geometry — constant for the lifetime of the finder.
  const s32 _loop_length;
  const s32 _sid_5p_length;
  const s32 _sid_3p_length;
  /// @brief Bitap search window for the loop: max(3×loop_length, 32) capped to the Bitap limit (64).
  const s32 _loop_search_window;

  // We can use Bitap to find the SID adapters as a fallback to the SIMD LUT-based search. As the bitap algorithm
  // finds the "end position" of the match, we will search in the forward direction for the 3p adapter. Only the 3p
  // bitaps are populated today; the matching 5p fallback path is not implemented.
  std::vector<Bitap<4>> _sid_3p_bitap;
  // Backward-direction counterpart of `_sid_3p_bitap`, used by `FindStartEnd` to recover the SID 3p
  // start position when falling back to the per-SID bitap search inside `CalculateMinimumErrorHairpin`.
  std::vector<Bitap<4>> _sid_3p_bw;

  // Helper struct to pass around some results of the duplex trim.
  struct TrimResults {
    u32 length{0};
  };

  // First step of duplex trim: look for the loop sequence using a string search and look for SIDs next to it.
  void FindHairpinByStringSearch(FixedReadRecord& record) const;

  // Second step of the duplex trim: after finding global symmetry in the read (very fast), look for the loop
  // and SID adapters around it.
  void FindHairpinByGlobalSymmetry(FixedReadRecord& record) const;

  // Third step: exhaustive search for the loop via local symmetry, then look for SID adapters if found.
  void FindHairpinByLocalSymmetry(FixedReadRecord& record, const TrimResults& results) const;

  /**
   * @brief Shared tail of every hairpin-search plan: run the SIMD SID search anchored to a known loop span.
   *
   * Uses two complementary anchor strategies in a single SIMD pass (4 of 8 lanes):
   * - Lanes 0–1 ("indel-aware"): anchored from loop_start, shifted by kSidMaxEditDistance to absorb
   *   SID drift from loop indels.
   * - Lanes 2–3 ("nominal"): anchored from loop_end using the nominal loop length, covering the
   *   majority case where the loop has no indels.
   *
   * The two strategies provide complementary coverage because the SIMD gather quantizes offsets to
   * 4-aligned bases; even a 1-base shift between them can move the quantized window by 4.
   *
   * @param[in,out] record Read record to populate with SID match / filter results.
   * @param[in] loop_start Start position (inclusive) of the matched loop k-mer in the read.
   * @param[in] loop_end End position (inclusive) of the matched loop k-mer in the read.
   */
  void FindSidsAroundLoop(FixedReadRecord& record, s32 loop_start, s32 loop_end) const;
};
}  // namespace xoos::demux
