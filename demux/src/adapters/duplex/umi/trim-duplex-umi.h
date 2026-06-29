#pragma once

#include "adapters/duplex/stem/trim-duplex-stem.h"
#include "io/read-record.h"

namespace xoos::demux {

class TrimDuplexUmi : public TrimDuplexStem {
 public:
  TrimDuplexUmi(SeqMatcher sid_5p_matcher, SeqMatcher sid_3p_matcher, SeqMatcher start_matcher, SeqMatcher end_matcher,
                const SeqMatcher& stem_5p_matcher, const SeqMatcher& stem_3p_matcher, SeqMatcher umi_5p_matcher,
                SeqMatcher umi_3p_matcher, std::string_view loop_sequence);

  std::pair<s32, s32> FindUMI(FixedReadRecord& record) const;

 private:
  SeqMatcher _umi_5p_matcher;
  SeqMatcher _umi_3p_matcher;

  // cascaded LUTs for UMI detection
  const CascadedLUTs _umi_lut;

  const s64 _mask_umi5p;
  const s64 _mask_umi3p;

  // These are for the SIMD LUT search, which searches up to 8 positions at once
  const s64 _lut_offset_umi_masks[8];
  const s64 _lut_offset_umi_types[8];

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
};

}  // namespace xoos::demux
