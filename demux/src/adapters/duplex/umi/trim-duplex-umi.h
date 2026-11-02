#pragma once

#include <array>

#include "adapters/duplex/duplex-match.h"
#include "adapters/duplex/stem/trim-duplex-stem.h"
#include "io/read-record.h"
#include "sequence/fixed-sequence.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

class TrimDuplexUmi : public TrimDuplexStem {
 public:
  TrimDuplexUmi(FixedSequence runway,
                FixedSequence start,
                const FixedSequence& stem,
                SeqMatcher umi_5p_matcher,
                SeqMatcher umi_3p_matcher);

  std::pair<s32, s32> FindUMI(FixedReadRecord& record) const;

 private:
  // These are for the SIMD LUT search, which searches up to 8 positions at once
  static constexpr size_t kSimdSearchLanes = 8;

  // The eight SIMD search offsets used to seed the UMI LUT search, derived from the bitap-anchored
  // stem positions (pos_5p/pos_3p) or theoretical positions when the stems are absent.
  std::array<int64_t, kSimdSearchLanes> ComputeSearchOffsets(const FixedReadRecord& record,
                                                             s32 pos_5p,
                                                             s32 pos_3p) const;

  SeqMatcher _umi_5p_matcher;
  SeqMatcher _umi_3p_matcher;

  // cascaded LUTs for UMI detection
  const CascadedLUTs _umi_lut;

  const s64 _mask_umi5p;
  const s64 _mask_umi3p;

  const std::array<s64, kSimdSearchLanes> _lut_offset_umi_masks;
  // Barcode type per SIMD lane; all values are compile-time constant.
  static constexpr std::array<s64, kSimdSearchLanes> _lut_offset_umi_types{DuplexMatch::BarcodeType::kUMI5p,
                                                                           DuplexMatch::BarcodeType::kUMI5p,
                                                                           DuplexMatch::BarcodeType::kUMI3p,
                                                                           DuplexMatch::BarcodeType::kUMI3p,
                                                                           DuplexMatch::BarcodeType::kUnknown,
                                                                           DuplexMatch::BarcodeType::kUnknown,
                                                                           DuplexMatch::BarcodeType::kUnknown,
                                                                           DuplexMatch::BarcodeType::kUnknown};

  // UMI lengths of actual region, without partial stem
  static constexpr s32 kUmiLen1 = 4;
  static constexpr s32 kUmiLen2 = 6;

  // Information about UMI lengths and offsets for partial matching
  struct UmiSearchInfo {
    // Lengths of the 5' (including part that is part of the stem to increase specificity of match)
    const s32 length_5p;

    const s32 umi_5p_offset_4bp{length_5p - kUmiLen1};
    const s32 umi_5p_offset_6bp{length_5p - kUmiLen2};
  };

  UmiSearchInfo _umi_info{UmiSearchInfo(static_cast<s32>(_umi_5p_matcher.Pool().front().sequence.length()))};

  const s64 _assumed_umi_len_5p_4bp{static_cast<s32>(RunwaySequence().length() + Start().length()) -
                                    _umi_info.umi_5p_offset_4bp};
  const s64 _assumed_umi_len_5p_6bp{static_cast<s32>(RunwaySequence().length() + Start().length()) -
                                    _umi_info.umi_5p_offset_6bp};
};

}  // namespace xoos::demux
