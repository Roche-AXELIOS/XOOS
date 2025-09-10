#include "xoos/yc-decode/yc-decoder.h"

#include <numeric>
#include <string>
#include <vector>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/util/string-functions.h>

namespace xoos::yc_decode {

// There are two decoding for discordant bases, one for each constituent read (R1 and R2) of the consensus sequence.
// Each character in the decoding is either a nucleotide or a gap ('-').
// Codes 'A', 'C', 'G', 'T' represent concordant bases. They are decoded to themselves.
//            discordant codes: ABCDEFGHIJKLMNOPQRSTUVWXYZ
constexpr char kDecodingR1[] = "ACCGTTGTAGG-A--C-ACT-GATC-";
constexpr char kDecodingR2[] = "AACACGGA--TAC-G-CGGT-CT-TT";
constexpr char kComplement[] = "TKGYRMCWXPBZF-QJOEVA-SHIDL";  // complement codes
const std::string kConcordantCodes = "ACGT";
const std::string kTagDelimiters = "+-";
constexpr char kTagDelimiterR1 = '+';
constexpr char kTagDelimiterR2 = '-';
constexpr char kDiscordantGapBase = '-';
const std::string kNumbers = "0123456789";

/**
 * @brief Deserialize YC tag from input string.
 * @param input_str input string containing YC tag, e.g. "123-45M67-89"
 * @return YC tag object
 */
YcTag DeserializeYcTag(const std::string& input_str) {
  // The YC tag string has 3 components delimited by either `+` or `-`:
  // 1. length of the left overhang
  // 2. duplex segment (for both R1 and R2), which may contain character codes for discordant bases
  // 3. length of the right overhang
  // Each overhang is assigned to either R1 or R2 depending on the associated delimiter.
  // An overhang with '+' belongs to R1, while an overhang with `-` belongs to R2.
  // Example:
  //   YC tag string:               "5+6E7-8"
  //   consensus sequence:   AAAAACCCCCCTCCCCCCCGGGGGGGG
  //   consensus quality:    33333~~~~~~!~~~~~~~33333333
  //   left overhang:        AAAAA
  //   duplex segment:            CCCCCCTCCCCCCC
  //   right overhang:                          GGGGGGGG
  //                                    v
  //   Decoded R1 sequence:  AAAAACCCCCCTCCCCCCC
  //   Decoded R2 sequence:       CCCCCCCCCCCCCCGGGGGGGG
  //                                    ^ code `E` in YC tag decodes to `T` in R1 and `C` in R2

  // Split the input string into components using the delimiters '+' and '-'
  const auto delimiter_pos1 = input_str.find_first_of(kTagDelimiters);
  const auto delimiter_pos2 = input_str.find_last_of(kTagDelimiters);

  if (delimiter_pos1 == std::string::npos || delimiter_pos2 == std::string::npos ||     // delimiter not found
      delimiter_pos1 == delimiter_pos2 ||                                               // delimiter found once
      input_str.find_first_of(kTagDelimiters, delimiter_pos1 + 1) != delimiter_pos2) {  // delimiter found >= 3x
    throw error::Error("YC tag '{}' is not in the correct format; expected 3 components delimited by `+` or `-`",
                       input_str);
  }
  const std::string left_overhang = input_str.substr(0, delimiter_pos1);
  const std::string duplex_segment = input_str.substr(delimiter_pos1 + 1, delimiter_pos2 - delimiter_pos1 - 1);
  const std::string right_overhang = input_str.substr(delimiter_pos2 + 1);
  u64 left_overhang_len = 0;
  u64 right_overhang_len = 0;
  // if left or right overhang is empty, it is assigned a length of 0
  // otherwise, it must be a number and must not contain any other characters
  if (!left_overhang.empty()) {
    if (left_overhang.find_first_not_of(kNumbers) != std::string::npos) {
      throw error::Error("YC tag '{}' is not in the correct format; left overhang length '{}' is not a number.",
                         input_str,
                         left_overhang);
    }
    left_overhang_len = std::stoul(left_overhang);
  }
  if (!right_overhang.empty()) {
    if (right_overhang.find_first_not_of(kNumbers) != std::string::npos) {
      throw error::Error("YC tag '{}' is not in the correct format; right overhang length '{}' is not a number.",
                         input_str,
                         right_overhang);
    }
    right_overhang_len = std::stoul(right_overhang);
  }
  vec<YcTagDuplexOp> duplex_ops{};
  u64 concordant_len = 0;
  for (const char c : duplex_segment) {
    if (c >= '0' && c <= '9') {
      // `c` is a digit
      if (concordant_len == 0) {
        concordant_len = c - '0';
      } else {
        // accumulate the number
        concordant_len = concordant_len * 10 + (c - '0');
      }
    } else if (c >= 'A' && c <= 'Z' && c != 'N' && c != 'U') {
      if (concordant_len > 0) {
        // store the previous concordant op
        duplex_ops.emplace_back(true, 0, concordant_len);
        concordant_len = 0;
      }
      // codes `ACGT` are concordant bases
      duplex_ops.emplace_back(kConcordantCodes.find(c) != std::string::npos, c, 1);
    } else {
      throw error::Error("Invalid code '{}' in YC tag '{}'", c, input_str);
    }
  }
  if (concordant_len > 0) {
    // store the final concordant op
    duplex_ops.emplace_back(true, 0, concordant_len);
  }
  return YcTag{left_overhang_len,
               duplex_ops,
               right_overhang_len,
               input_str[delimiter_pos1] == kTagDelimiterR1,
               input_str[delimiter_pos2] == kTagDelimiterR1};
}

/**
 * @brief Deserialize YC tag from BAM record, reverse-complement and trim it as necessary.
 * @param record BAM record pointer
 * @return YC tag object
 * @details If the sequence is aligned to the reverse strand, then the YC tag would be reverse-complemented.
 * YC tag that is longer than the sequence may be trimmed based on left and right hard clips.
 */
YcTag DeserializeYcTag(const bam1_t* record) {
  // extract YC tag string from BAM record aux fields
  const u8* yc_tag_data = bam_aux_get(record, "YC");
  if (yc_tag_data == nullptr) {
    return YcTag{};  // return empty YC tag
  }
  const auto yc_tag_str = std::string(bam_aux2Z(yc_tag_data));
  // extract hard clip lengths from CIGAR string
  u64 left_hard_clip = 0;
  u64 right_hard_clip = 0;
  if (const auto cigar_ops = record->core.n_cigar; cigar_ops > 1) {
    const u32* cigar = bam_get_cigar(record);
    if (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP) {
      left_hard_clip = bam_cigar_oplen(cigar[0]);
    }
    if (bam_cigar_op(cigar[cigar_ops - 1]) == BAM_CHARD_CLIP) {
      right_hard_clip = bam_cigar_oplen(cigar[cigar_ops - 1]);
    }
  }
  // deserialize YC tag string and process it as necessary
  auto yc_tag = DeserializeYcTag(yc_tag_str);
  if (bam_is_rev(record)) {
    // if reverse-complement is required, it must be performed before trimming
    yc_tag.ReverseComplement();
  }
  const auto seq_len = static_cast<u64>(record->core.l_qseq);
  if ((left_hard_clip > 0 || right_hard_clip > 0) && seq_len < yc_tag.GetSequenceLength()) {
    // BAM record has hard clip(s), and the YC tag is longer than the consensus sequence.
    // Trim the YC tag based on the hard clips.
    yc_tag.Trim(left_hard_clip, right_hard_clip);
  }
  return yc_tag;
}

/**
 * @brief Reverse-complement YC tag in place.
 */
void YcTag::ReverseComplement() {
  // swap left and right overhangs and their read assignments, but do not flip the read assignments (e.g. R1 -> R2)
  // because reverse-complemented R1 and R2 are still R1 and R2
  std::swap(left_overhang, right_overhang);
  std::swap(left_overhang_is_r1, right_overhang_is_r1);
  if (!duplex_ops.empty()) {
    // reverse the duplex ops and complement any valid codes
    std::ranges::reverse(duplex_ops);
    for (auto& op : duplex_ops) {
      if (op.code >= 'A' && op.code <= 'Z' && op.code != 'N' && op.code != 'U') {
        op.code = kComplement[op.code - 'A'];
      }
    }
  }
}

/**
 * @brief Return the total sequence length of the YC tag.
 * @return Total sequence length
 */
u64 YcTag::GetSequenceLength() const {
  u64 length = left_overhang + right_overhang;
  for (const auto& op : duplex_ops) {
    length += op.length;
  }
  return length;
}

/**
 * @brief Trim both ends of YC tag in place.
 * @param left_clip Left clip length
 * @param right_clip Right clip length
 */
void YcTag::Trim(u64 left_clip, u64 right_clip) {
  const u64 seq_len = this->GetSequenceLength();
  const u64 duplex_segment_length = seq_len - left_overhang - right_overhang;
  const u64 target_len = seq_len - left_clip - right_clip;
  // process the left overhang
  if (left_clip < left_overhang) {
    // trim overhang and consume entire clip length
    left_overhang -= left_clip;
    left_clip = 0;
  } else {
    // consume entire overhang and reduce clip length
    left_clip -= left_overhang;
    left_overhang = 0;
  }
  // process the right overhang
  if (right_clip < right_overhang) {
    right_overhang -= right_clip;
    right_clip = 0;
  } else {
    right_clip -= right_overhang;
    right_overhang = 0;
  }
  // process the duplex segment
  if (!duplex_ops.empty() && (left_clip > 0 || right_clip > 0)) {
    vec<YcTagDuplexOp> new_duplex_ops{};
    if (left_clip >= duplex_segment_length) {
      // left clip trims the entire duplex segment
      left_clip -= duplex_segment_length;
    } else if (right_clip >= duplex_segment_length) {
      // right clip trims the entire duplex segment
      right_clip -= duplex_segment_length;
    } else {
      // left/right clips partially trim the duplex segment
      u64 idx = left_overhang;
      for (const auto& [concordant, code, length] : duplex_ops) {
        auto op_len = length;
        if (left_clip > 0) {
          if (left_clip < op_len) {
            // trim op by left clip length
            op_len -= left_clip;
            left_clip = 0;
          } else {
            // reduce left clip length
            left_clip -= op_len;
            op_len = 0;
          }
        }
        if (left_clip == 0) {
          const bool has_right_clip = idx + op_len > target_len;
          if (has_right_clip) {
            // reduce op length because it is partially trimmed by right clip
            op_len = target_len - idx;
          }
          if (op_len > 0) {
            idx += op_len;
            new_duplex_ops.emplace_back(concordant, code, op_len);
          }
          if (has_right_clip) {
            // stop here because we are at the op partially trimmed by right clip
            // the rest of the duplex segment is completely trimmed by right clip
            break;
          }
        }
      }
    }
    duplex_ops = new_duplex_ops;
  }
  if (duplex_ops.empty()) {
    // the entire duplex segment was trimmed
    if (left_clip > 0) {
      // trim right overhang with left clip
      if (right_overhang > left_clip) {
        right_overhang -= left_clip;
      } else {
        right_overhang = 0;
      }
    }
    if (right_clip > 0) {
      // trim left overhang with right clip
      if (left_overhang > right_clip) {
        left_overhang -= right_clip;
      } else {
        left_overhang = 0;
      }
    }
  }
}

/**
 * @brief Check whether the consensus read is duplex.
 * @return True if the read is duplex.
 */
bool YcTag::IsDuplex() const {
  return !duplex_ops.empty();
}

/**
 * @brief Check whether the consensus read is simplex. Special case: an empty YC tag is considered simplex.
 * @return True if the read is simplex.
 */
bool YcTag::IsSimplex() const {
  return duplex_ops.empty();
}

/**
 * @brief Count the number of discordant bases in the duplex segment.
 * @return Number of discordant bases
 */
u64 YcTag::CountDiscordantBase() const {
  // Note that "simplex" reads have no discordant bases
  return std::ranges::count_if(duplex_ops, [](const YcTagDuplexOp& op) { return !op.is_numeric; });
}

/**
 * @brief Count the number of bases in the duplex segment, including both concordant and discordant bases.
 * @return Number of duplex bases
 */
u64 YcTag::CountDuplexBase() const {
  return std::accumulate(
      duplex_ops.begin(), duplex_ops.end(), 0UL, [](u64 sum, const YcTagDuplexOp& op) { return sum + op.length; });
}

/**
 * @brief Extract the base type at each position represented by the YC tag.
 * @return Vector of base type enums
 */
vec<BaseType> YcTag::GetBaseTypes() const {
  using enum BaseType;
  vec<BaseType> base_categories{};
  base_categories.reserve(GetSequenceLength());
  if (left_overhang > 0) {
    base_categories.insert(base_categories.begin(), left_overhang, kSimplex);
  }
  for (const auto& op : duplex_ops) {
    if (op.is_numeric) {
      base_categories.insert(base_categories.end(), op.length, kConcordant);
    } else {
      base_categories.insert(base_categories.end(), op.length, kDiscordant);
    }
  }
  if (right_overhang > 0) {
    base_categories.insert(base_categories.end(), right_overhang, kSimplex);
  }
  return base_categories;
}

/**
 * @brief Serialize the YC tag to string format.
 * @return string format of the YC tag
 */
std::string YcTag::ToString() const {
  std::string yc_tag_str{};
  yc_tag_str += std::to_string(left_overhang);
  if (left_overhang_is_r1) {
    yc_tag_str += kTagDelimiterR1;
  } else {
    yc_tag_str += kTagDelimiterR2;
  }
  for (const auto& op : duplex_ops) {
    if (op.is_numeric) {
      yc_tag_str += std::to_string(op.length);
    } else {
      yc_tag_str += op.code;
    }
  }
  if (right_overhang_is_r1) {
    yc_tag_str += kTagDelimiterR1;
  } else {
    yc_tag_str += kTagDelimiterR2;
  }
  yc_tag_str += std::to_string(right_overhang);
  return yc_tag_str;
}

/**
 * @brief Encode base types to YC tag string format.
 * @param base_types Base types
 * @param reverse Flag whether to reverse the order of base types
 * @return YC tag string
 */
std::string EncodeBaseTypesToYcTag(const vec<BaseType>& base_types, bool reverse) {
  using enum BaseType;
  u64 left_overhang = 0;
  std::string duplex_segment{};
  u64 num_simplex = 0;
  u64 num_concordant = 0;
  auto process_type = [&](const BaseType& base_type) {
    switch (base_type) {
      case kDiscordant:
        if (num_simplex > 0) {
          left_overhang = num_simplex;
          num_simplex = 0;
        }
        if (num_concordant > 0) {
          duplex_segment.append(std::to_string(num_concordant));
          num_concordant = 0;
        }
        duplex_segment += 'X';  // 'X' was chosen arbitrarily; any discordant character code would also work.
        break;
      case kConcordant:
        if (num_simplex > 0) {
          left_overhang = num_simplex;
          num_simplex = 0;
        }
        ++num_concordant;
        break;
      case kSimplex:
        ++num_simplex;
        break;
    }
  };
  if (reverse) {
    // YC tag needs to be reverse complemented
    std::for_each(base_types.rbegin(), base_types.rend(), process_type);
  } else {
    std::ranges::for_each(base_types, process_type);
  }
  if (num_concordant > 0) {
    // append the last concordant op
    duplex_segment.append(std::to_string(num_concordant));
  }
  if (left_overhang == 0 && duplex_segment.empty() && num_simplex > 0) {
    // no duplex segment, but simplex bases are present
    left_overhang = num_simplex;
    num_simplex = 0;
  }
  return fmt::format("{}+{}+{}", left_overhang, duplex_segment, num_simplex);
}

/**
 * @brief Decode the YC tag string for a consensus sequence into two constituent reads.
 * Consensus sequence, base quality, and YC tag are expected to be in the same orientation.
 * Decoded reads will be in the same orientation as the consensus sequence.
 * @param consensus_seq Consensus sequence
 * @param consensus_qual Consensus quality
 * @param yc_tag YC tag
 * @return Decoded reads
 */
DecodedReads Decode(const std::string& consensus_seq, const std::string& consensus_qual, const YcTag& yc_tag) {
  const auto seq_len = consensus_seq.length();
  if (seq_len != consensus_qual.length()) {
    throw error::Error(
        "Consensus sequence '{}' and quality '{}' have unequal lengths", seq_len, consensus_qual.length());
  }
  const auto yc_tag_len = yc_tag.GetSequenceLength();
  if (seq_len != yc_tag_len) {
    throw error::Error("Consensus sequence '{}' and YC tag '{}' have unequal lengths", seq_len, yc_tag_len);
  }
  // Both output reads would be in the same orientation as the consensus sequence.
  std::string r1_seq{};
  std::string r2_seq{};
  std::string r1_qual{};
  std::string r2_qual{};
  using enum BaseType;
  vec<BaseType> r1_base_types{};
  vec<BaseType> r2_base_types{};
  // process the left overhang
  if (yc_tag.left_overhang > 0) {
    if (yc_tag.left_overhang_is_r1) {
      // append the simplex segment to R1
      r1_seq = consensus_seq.substr(0, yc_tag.left_overhang);
      r1_qual = consensus_qual.substr(0, yc_tag.left_overhang);
      r1_base_types.insert(r1_base_types.end(), yc_tag.left_overhang, kSimplex);
    } else {
      // append the simplex segment to R2
      r2_seq = consensus_seq.substr(0, yc_tag.left_overhang);
      r2_qual = consensus_qual.substr(0, yc_tag.left_overhang);
      r2_base_types.insert(r2_base_types.end(), yc_tag.left_overhang, kSimplex);
    }
  }
  // process the duplex segment
  u64 idx = yc_tag.left_overhang;
  for (const auto& [concordant, code, length] : yc_tag.duplex_ops) {
    if (concordant) {
      if (length > 0) {
        const auto seq = consensus_seq.substr(idx, length);
        r1_seq += seq;
        r2_seq += seq;
        const auto qual = consensus_qual.substr(idx, length);
        r1_qual += qual;
        r2_qual += qual;
        r1_base_types.insert(r1_base_types.end(), length, kConcordant);
        r2_base_types.insert(r2_base_types.end(), length, kConcordant);
        idx += length;
      }
    } else {
      // this is a discordant base
      const auto qual = consensus_qual.at(idx);
      const auto r1_base = kDecodingR1[code - 'A'];
      const auto r2_base = kDecodingR2[code - 'A'];
      if (r1_base == kDiscordantGapBase) {
        if (!r1_base_types.empty()) {
          // The gap does not add a base (or quality score or base type) to R1, but a base is added to R2 instead.
          // To indicate this for R1, modify the base type at the previous position in R1 to "discordant".
          // A distinction between discordant and concordant indels is crucial to accurate variant calling.
          r1_base_types.back() = kDiscordant;
        }
      } else {
        r1_seq += r1_base;
        r1_qual += qual;
        r1_base_types.insert(r1_base_types.end(), 1, kDiscordant);
      }
      if (r2_base == kDiscordantGapBase) {
        if (!r2_base_types.empty()) {
          // Modify the base type at the previous position to "discordant".
          r2_base_types.back() = kDiscordant;
        }
      } else {
        r2_seq += r2_base;
        r2_qual += qual;
        r2_base_types.insert(r2_base_types.end(), 1, kDiscordant);
      }
      ++idx;
    }
  }
  if (idx < seq_len) {
    if (yc_tag.right_overhang_is_r1) {
      // append the remaining simplex bases to R1
      r1_seq += consensus_seq.substr(idx);
      r1_qual += consensus_qual.substr(idx);
      r1_base_types.insert(r1_base_types.end(), yc_tag.right_overhang, kSimplex);
    } else {
      // append the remaining simplex bases to R2
      r2_seq += consensus_seq.substr(idx);
      r2_qual += consensus_qual.substr(idx);
      r2_base_types.insert(r2_base_types.end(), yc_tag.right_overhang, kSimplex);
    }
  }
  return DecodedReads{r1_seq, r2_seq, r1_qual, r2_qual, r1_base_types, r2_base_types};
}

/**
 * @brief Decode the YC tag string for a consensus sequence into two constituent reads.
 * Consensus sequence and base quality are expected to be in the same orientation.
 * YC tag in the opposite orientation will be reverse-complemented.
 * YC tag that is longer than the consensus sequence may be trimmed based on left/right clips.
 * Decoded reads will be in the same orientation as the consensus sequence.
 * @param consensus_seq Consensus sequence
 * @param consensus_qual Consensus quality
 * @param yc_tag_str YC tag string
 * @param rev_comp Flag whether YC tag needs to be reverse-complemented
 * @param left_clip Left clip length for trimming YC tag
 * @param right_clip Right clip length for trimming YC tag
 * @return Decoded reads
 */
DecodedReads Decode(const std::string& consensus_seq,
                    const std::string& consensus_qual,
                    const std::string& yc_tag_str,
                    const bool rev_comp,
                    const u64 left_clip,
                    const u64 right_clip) {
  YcTag yc_tag = DeserializeYcTag(yc_tag_str);
  if (rev_comp) {
    // If reverse-complement is required, it must be done before trimming
    yc_tag.ReverseComplement();
  }
  if ((left_clip > 0 || right_clip > 0) && consensus_seq.size() < yc_tag.GetSequenceLength()) {
    // There are hard clips, and the YC tag is longer than the consensus sequence.
    // Trim the YC tag based on the hard clips.
    yc_tag.Trim(left_clip, right_clip);
  }
  return Decode(consensus_seq, consensus_qual, yc_tag);
}

/**
 * @brief Helper function to extract alignment information for decoded alignments.
 * @param yc_tag YC tag
 * @param record BAM record
 * @param decoded_alns Decoded alignments to be updated with alignment information
 */
static void UpdateAlignmentInfo(const YcTag& yc_tag, const bam1_t* record, DecodedAlignments& decoded_alns) {
  if (auto cigar_infos = GetCigarAlignInfos(record); !cigar_infos.empty()) {
    u64 cigar_info_idx = 0;
    u64 yc_tag_op_read_start = 0;

    const bool has_duplex = !yc_tag.duplex_ops.empty();
    u64 r1_ref_pos = 0;
    u64 r2_ref_pos = 0;
    if (yc_tag.left_overhang_is_r1) {
      r1_ref_pos = record->core.pos;
      if (has_duplex) {
        r2_ref_pos = r1_ref_pos;
      }
    } else {
      r2_ref_pos = record->core.pos;
      if (has_duplex) {
        r1_ref_pos = r2_ref_pos;
      }
    }
    vec<HtslibCigarOp> r1_cigar;
    vec<HtslibCigarOp> r2_cigar;
    // process CIGAR op(s) within the left overhang
    if (yc_tag.left_overhang > 0) {
      const auto simplex_tail = yc_tag.left_overhang;
      for (; cigar_info_idx < cigar_infos.size(); ++cigar_info_idx) {
        auto [cigar_op, read_pos, ref_pos, read_span, ref_span] = cigar_infos.at(cigar_info_idx);
        if (read_pos >= simplex_tail) {
          // CIGAR op starts after the overhang
          break;
        }
        auto& r_cigar = yc_tag.left_overhang_is_r1 ? r1_cigar : r2_cigar;
        if (read_pos + read_span > simplex_tail) {
          // CIGAR op spans the end of the overhang and beyond
          const u64 span = simplex_tail - read_pos;
          r_cigar.emplace_back(cigar_op, span);
          if (has_duplex) {
            if (yc_tag.left_overhang_is_r1) {
              r2_ref_pos = ref_span > 0 ? ref_pos + span : ref_pos;
            } else {
              r1_ref_pos = ref_span > 0 ? ref_pos + span : ref_pos;
            }
          }
          break;
        }
        // CIGAR op is consumed completely by the overhang
        r_cigar.emplace_back(cigar_op, std::max(read_span, ref_span));
        if (has_duplex) {
          if (yc_tag.left_overhang_is_r1) {
            r2_ref_pos = ref_pos + ref_span;
          } else {
            r1_ref_pos = ref_pos + ref_span;
          }
        }
      }
      yc_tag_op_read_start = yc_tag.left_overhang;
    }
    // process CIGAR op(s) within the duplex segment
    if (has_duplex) {
      // if left overhang is empty, then the duplex segment starts at the beginning of the read
      // adjust reference positions if the leading YC tag ops decode to gap(s) in R1/R2 with respect to the consensus
      bool adjust_r1_ref_pos = !yc_tag.left_overhang_is_r1 || yc_tag.left_overhang == 0;
      bool adjust_r2_ref_pos = yc_tag.left_overhang_is_r1 || yc_tag.left_overhang == 0;
      for (const auto& yc_tag_op : yc_tag.duplex_ops) {
        const auto yc_tag_op_len = yc_tag_op.length;
        auto yc_tag_op_read_end = yc_tag_op_read_start + yc_tag_op_len;
        for (; cigar_info_idx < cigar_infos.size(); ++cigar_info_idx) {
          auto [cigar_op, read_pos, ref_pos, read_span, ref_span] = cigar_infos.at(cigar_info_idx);
          if (read_pos >= yc_tag_op_read_end) {
            // CIGAR op starts after the duplex segment
            break;
          }
          const bool cigar_op_consumes_read = ConsumesRead(cigar_op);
          const bool cigar_op_consumes_ref = ConsumesReference(cigar_op);
          if (cigar_op_consumes_ref && !cigar_op_consumes_read) {
            // CIGAR op consumes reference base(s) but not read base(s)
            adjust_r1_ref_pos = false;
            adjust_r2_ref_pos = false;
            if (!r1_cigar.empty() && r1_cigar.back().op == cigar_op) {
              r1_cigar.back().len += ref_span;
            } else {
              r1_cigar.emplace_back(cigar_op, ref_span);
            }
            if (!r2_cigar.empty() && r2_cigar.back().op == cigar_op) {
              r2_cigar.back().len += ref_span;
            } else {
              r2_cigar.emplace_back(cigar_op, ref_span);
            }
          } else if (yc_tag_op.is_numeric) {
            adjust_r1_ref_pos = false;
            adjust_r2_ref_pos = false;
            // CIGAR op consumes read base(s)
            // Determine the span of CIGAR op with respect to the YC tag op's read start and end positions
            const u64 span =
                std::min(read_pos + read_span, yc_tag_op_read_end) - std::max(read_pos, yc_tag_op_read_start);
            if (!r1_cigar.empty() && r1_cigar.back().op == cigar_op) {
              r1_cigar.back().len += span;
            } else {
              r1_cigar.emplace_back(cigar_op, span);
            }
            if (!r2_cigar.empty() && r2_cigar.back().op == cigar_op) {
              r2_cigar.back().len += span;
            } else {
              r2_cigar.emplace_back(cigar_op, span);
            }
          } else {
            // CIGAR op consumes a read base
            // since YC tag op is discordant, decode the character code
            const auto decoding_idx = yc_tag_op.code - 'A';
            if (kDecodingR1[decoding_idx] == kDiscordantGapBase) {
              // R1 base is a gap with respect to the consensus sequence
              if (cigar_op_consumes_ref) {
                // since CIGAR op also consumes the reference, we must add a deletion CIGAR op to reflect the gap
                if (adjust_r1_ref_pos) {
                  ++r1_ref_pos;
                }
                if (!r1_cigar.empty() && r1_cigar.back().op == BAM_CDEL) {
                  r1_cigar.back().len += yc_tag_op_len;
                } else {
                  r1_cigar.emplace_back(BAM_CDEL, yc_tag_op_len);
                }
              }
            } else {
              // R1 base is not a gap base
              adjust_r1_ref_pos = false;
              if (!r1_cigar.empty() && r1_cigar.back().op == cigar_op) {
                r1_cigar.back().len += yc_tag_op_len;
              } else {
                r1_cigar.emplace_back(cigar_op, yc_tag_op_len);
              }
            }
            if (kDecodingR2[decoding_idx] == kDiscordantGapBase) {
              // R2 base is a gap with respect to the consensus sequence
              if (cigar_op_consumes_ref) {
                // since CIGAR op also consumes the reference, we must add a deletion CIGAR op to reflect the gap
                if (adjust_r2_ref_pos) {
                  ++r2_ref_pos;
                }
                if (!r2_cigar.empty() && r2_cigar.back().op == BAM_CDEL) {
                  r2_cigar.back().len += yc_tag_op_len;
                } else {
                  r2_cigar.emplace_back(BAM_CDEL, yc_tag_op_len);
                }
              }
            } else {
              // R2 base is not a gap base
              adjust_r2_ref_pos = false;
              if (!r2_cigar.empty() && r2_cigar.back().op == cigar_op) {
                r2_cigar.back().len += yc_tag_op_len;
              } else {
                r2_cigar.emplace_back(cigar_op, yc_tag_op_len);
              }
            }
          }
          if (read_pos + read_span > yc_tag_op_read_end) {
            // CIGAR op spans beyond the duplex segment, do not look at the next CIGAR op
            break;
          }
        }
        yc_tag_op_read_start = yc_tag_op_read_end;
      }
    }
    // process CIGAR op(s) within the right overhang
    if (yc_tag.right_overhang > 0) {
      auto [cigar_op, read_pos, ref_pos, read_span, ref_span] = cigar_infos.at(cigar_info_idx);
      auto& r_cigar = yc_tag.right_overhang_is_r1 ? r1_cigar : r2_cigar;
      if (read_pos < yc_tag_op_read_start) {
        // CIGAR op was partially consumed in the duplex segment, append the remaining part of CIGAR op
        const u64 span = read_pos + read_span - yc_tag_op_read_start;
        if (!r_cigar.empty() && r_cigar.back().op == cigar_op) {
          r_cigar.back().len += span;
        } else {
          r_cigar.emplace_back(cigar_op, span);
        }
        ++cigar_info_idx;
      }
      // append all remaining CIGAR ops
      for (; cigar_info_idx < cigar_infos.size(); ++cigar_info_idx) {
        const auto& info = cigar_infos.at(cigar_info_idx);
        r_cigar.emplace_back(info.op, std::max(info.read_span, info.ref_span));
      }
    }
    // remove leading/trailing deletions and convert leading/trailing insertions to soft-clips
    if (!r1_cigar.empty()) {
      const auto& last_op = r1_cigar.back().op;
      if (last_op == BAM_CDEL) {
        r1_cigar.pop_back();
      } else if (last_op == BAM_CINS) {
        r1_cigar.back().op = BAM_CSOFT_CLIP;
      }
    }
    if (!r2_cigar.empty()) {
      const auto& last_op = r2_cigar.back().op;
      if (last_op == BAM_CDEL) {
        r2_cigar.pop_back();
      } else if (last_op == BAM_CINS) {
        r2_cigar.back().op = BAM_CSOFT_CLIP;
      }
    }
    if (!r1_cigar.empty()) {
      const auto& first_op = r1_cigar.front().op;
      if (first_op == BAM_CDEL) {
        r1_cigar.erase(r1_cigar.begin());
      } else if (first_op == BAM_CINS) {
        r1_cigar.front().op = BAM_CSOFT_CLIP;
      }
    }
    if (!r2_cigar.empty()) {
      const auto& first_op = r2_cigar.front().op;
      if (first_op == BAM_CDEL) {
        r2_cigar.erase(r2_cigar.begin());
      } else if (first_op == BAM_CINS) {
        r2_cigar.front().op = BAM_CSOFT_CLIP;
      }
    }
    decoded_alns.r1_ref_pos = r1_ref_pos;
    decoded_alns.r2_ref_pos = r2_ref_pos;
    decoded_alns.r1_cigar = r1_cigar;
    decoded_alns.r2_cigar = r2_cigar;
  }
}

/**
 * @brief Decode a BAM record using its YC tag.
 * Decoded alignments will be in the same orientation as the consensus sequence from the BAM record.
 * @param record BAM record
 * @return Decoded alignments
 */
DecodedAlignments Decode(const bam1_t* record) {
  auto yc_tag = DeserializeYcTag(record);
  if (yc_tag.GetSequenceLength() == 0) {
    return DecodedAlignments{};
  }
  // extract consensus sequence and base quality
  std::string consensus_seq;
  std::string consensus_qual;
  const auto seq_len = record->core.l_qseq;
  consensus_seq.resize(seq_len);
  consensus_qual.resize(seq_len);
  auto* const seq_data = bam_get_seq(record);
  auto* const qual_data = bam_get_qual(record);
  for (int i = 0; i < seq_len; ++i) {
    consensus_seq[i] = seq_nt16_str[bam_seqi(seq_data, i)];
    consensus_qual[i] = static_cast<char>(qual_data[i] + kPhred33Offset);
  }
  // decode the consensus sequence using the YC tag
  auto [r1_seq, r2_seq, r1_qual, r2_qual, r1_base_types, r2_base_types] = Decode(consensus_seq, consensus_qual, yc_tag);
  DecodedAlignments decoded_alns;
  decoded_alns.r1_seq = r1_seq;
  decoded_alns.r2_seq = r2_seq;
  decoded_alns.r1_qual = r1_qual;
  decoded_alns.r2_qual = r2_qual;
  decoded_alns.r1_base_types = r1_base_types;
  decoded_alns.r2_base_types = r2_base_types;
  UpdateAlignmentInfo(yc_tag, record, decoded_alns);
  return decoded_alns;
}

/**
 * @brief Add suffix to the read name. The suffix is added to the first string before `:` in the name.
 * @param name Read name
 * @param suffix Suffix to be added
 */
static std::string AddNameSuffixHelper(const std::string& name, const std::string& suffix) {
  if (const u64 pos = name.find(':'); pos != std::string::npos) {
    return name.substr(0, pos) + suffix + name.substr(pos);
  }
  return name + suffix;
}

// Helper function to set up a new BAM record based on the source BAM record
static void DecodeToBamRecordsHelper(const bam1_t* src_record,
                                     const io::Bam1Ptr& new_record,
                                     const std::string& seq,
                                     const std::string& qual,
                                     const hts_pos_t ref_pos,
                                     const std::string& qname,
                                     const vec<HtslibCigarOp>& cigar_ops,
                                     const std::string& yc_tag) {
  // set up the cigar string data
  vec<u32> cigar;
  for (const auto& [op, len] : cigar_ops) {
    cigar.push_back(bam_cigar_gen(len, op));
  }
  // shift quality values to Phred33
  std::string qual_shifted;
  qual_shifted.resize(qual.size());
  for (u64 i = 0; i < qual.size(); ++i) {
    qual_shifted[i] = static_cast<char>(qual[i] - kPhred33Offset);
  }
  // set up the new BAM record
  const auto aux_data_len = bam_get_l_aux(src_record);
  bam_set1(new_record.get(),
           qname.size(),
           qname.c_str(),
           src_record->core.flag,
           src_record->core.tid,
           ref_pos,
           src_record->core.qual,
           cigar.size(),
           cigar.data(),
           src_record->core.mtid,
           src_record->core.mpos,
           src_record->core.isize,
           seq.size(),
           seq.data(),
           qual_shifted.data(),
           aux_data_len);
  // Iterate over the auxiliary fields (BAM tags) of the source record and naively copy them to the new record.
  // If the consensus sequence contains discordant bases (character codes in YC tag), then the auxiliary fields
  // (e.g. `NM` and `MD`) are very likely to be incorrect for the decoded reads. For example, a discordant code may
  // convert a reference mismatch (`NM:i:1`) from the consensus sequence into a reference match (`NM:i:0`) in a decoded
  // read.
  // TODO: adjust the auxiliary data to the correct values
  const u8* aux_first = bam_aux_first(src_record);
  for (const u8* aux = aux_first; aux != nullptr;) {
    const u8* next_aux = bam_aux_next(src_record, aux);
    // Determine the length of the auxiliary data
    const s64 length = (next_aux != nullptr) ? next_aux - aux : aux_first + aux_data_len - aux;
    // Append the auxiliary data to the destination record
    // Subtract 3 from the length to account for the tag (2 bytes) and the type (1 byte)
    bam_aux_append(new_record.get(), bam_aux_tag(aux), bam_aux_type(aux), static_cast<int>(length) - 3, aux + 1);
    aux = next_aux;
  }
  bam_aux_update_str(new_record.get(), "YC", static_cast<s32>(yc_tag.length() + 1), yc_tag.c_str());
}

/**
 * @brief Decode a BAM record using its YC tag.
 * If the sequence is aligned to the reverse strand, then the YC tag would be reverse-complemented.
 * YC tag that is longer than the sequence may be trimmed based on left and right hard clips.
 * Decoded alignments will be in the same orientation as the consensus sequence from the BAM record.
 * Left and right overhangs, if present, would be assigned to R1.
 * Output BAM records are generated for R1 and R2 based on the input BAM record with the following adjusted:
 * 1. read name
 * 2. reference start position
 * 3. CIGAR string
 * 4. read sequence
 * 5. base quality
 * @param record BAM record
 * @return one of the following:
 * 1. empty, if no reads can be decoded
 * 2. one BAM record, if there is a overhang but no duplex segment
 * 3. two BAM records, if there is a duplex segment
 */
std::variant<std::monostate, io::Bam1Ptr, std::pair<io::Bam1Ptr, io::Bam1Ptr>> DecodeToBamRecords(
    const bam1_t* record) {
  const bool rev_comp = bam_is_rev(record);
  const auto decoded_alns = Decode(record);
  const bool has_r1 = !decoded_alns.r1_seq.empty();
  const bool has_r2 = !decoded_alns.r2_seq.empty();
  if (!has_r1 && !has_r2) {
    return std::monostate{};
  }
  const auto qname = std::string{bam_get_qname(record)};
  auto r1_record = io::Bam1Ptr(bam_init1());
  if (has_r1) {
    // Set up the BAM record for R1 using the original record
    // The YC tag is updated to encode for base types in R1
    const auto r1_qname = AddNameSuffixHelper(qname, "/1");
    DecodeToBamRecordsHelper(record,
                             r1_record,
                             decoded_alns.r1_seq,
                             decoded_alns.r1_qual,
                             static_cast<hts_pos_t>(decoded_alns.r1_ref_pos),
                             r1_qname,
                             decoded_alns.r1_cigar,
                             EncodeBaseTypesToYcTag(decoded_alns.r1_base_types, rev_comp));
    if (!has_r2) {
      return std::move(r1_record);
    }
  }
  auto r2_record = io::Bam1Ptr(bam_init1());
  if (has_r2) {
    // Set up the BAM record for R2 using the original record
    // The YC tag is updated to encode for base types in R2
    const auto r2_qname = AddNameSuffixHelper(qname, "/2");
    DecodeToBamRecordsHelper(record,
                             r2_record,
                             decoded_alns.r2_seq,
                             decoded_alns.r2_qual,
                             static_cast<hts_pos_t>(decoded_alns.r2_ref_pos),
                             r2_qname,
                             decoded_alns.r2_cigar,
                             EncodeBaseTypesToYcTag(decoded_alns.r2_base_types, rev_comp));
    if (!has_r1) {
      return std::move(r2_record);
    }
  }
  return std::pair{std::move(r1_record), std::move(r2_record)};
}

/**
 * Constructs an std::runtime_error indicating an invalid base type in the read with the provided read name.
 *
 * @param read_name The name of the read where the invalid base type was encountered.
 * @return std::runtime_error The constructed runtime error with the formatted message.
 */
std::runtime_error BaseTypeError(const std::string& read_name) {
  return error::Error("Invalid base type in read named '{}'", read_name);
}

}  // namespace xoos::yc_decode
