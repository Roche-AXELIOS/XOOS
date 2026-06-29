#include "xoos/yc-decode/yc-decoder.h"

#include <numeric>
#include <string>
#include <vector>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/util/string-functions.h>

namespace xoos::yc_decode {

// Each character code in the YC tag indicates a single differing base between the consensus sequence and R1/R2.
// Upper case codes represent mismatches and insertions in the consensus sequence with respect to R1 or R2. The
// consensus has a base at that position while either R1 or R2 has a different base or a gap.
// Lower case codes represent deletions in the consensus sequence with respect to R1 or R2. The consensus has a gap at
// that position while either R1 or R2 has a base.
// Codes 'U', 'N', 'u', and 'n` are not supported.
// The decoding contains either a nucleotide ('A', 'C', 'G', 'T') or a gap ('-').
//            upper case codes: ABCDEFGHIJKLMNOPQRSTUVWXYZ
//            lower case codes: abcdefghijklmnopqrstuvwxyz
constexpr char kDecodingR1[] = "ACCGTTGTAGG-A--C-ACT-GATC-";
constexpr char kDecodingR2[] = "AACACGGA--TAC-G-CGGT-CT-TT";
constexpr char kComplement[] = "TKGYRMCWXPBZF-QJOEVA-SHIDL";
constexpr char kGapComplement[] = "tkgyrmcwxpbzf-qjoeva-shidl";

const std::string kTagDelimiters = "+-";
constexpr char kTagDelimiterR1 = '+';
constexpr char kTagDelimiterR2 = '-';
constexpr char kDiscordantGapBase = '-';
const std::string kNumbers = "0123456789";

// Char for discordant base type encoding to YC tag.
// 'X' was chosen arbitrarily; any discordant character code would also work.
constexpr char kDiscordantChar = 'X';

// Gap base in the consensus creates an inserted base in R1/R2 with no known base quality value, set it to 0 (!).
const std::string kConsensusGapBaseQuality = "!";

/**
 * @brief Check if a code encodes a consensus mismatch to either R1 or R2.
 * @param code Discordant character code
 * @return True if the consensus base is a mismatch to R1 or R2, false otherwise.
 */
static bool IsMismatchToR1R2(const char code) {
  return code >= 'A' && code <= 'Z' && code != 'N' && code != 'U';
}

/**
 * @brief Check if a code encodes a consensus gap to either R1 or R2.
 * @param code Discordant character code
 * @return True if the consensus base is a gap to R1 or R2, false otherwise.
 */
static bool IsGapToR1R2(const char code) {
  return code >= 'a' && code <= 'z' && code != 'n' && code != 'u';
}

/**
 * @brief Decode a character code to a base in R1 or R2 using the provided decoding table.
 * @param code Discordant character code
 * @param decoding_table Decoding table for R1 or R2
 * @return Decoded base in R1 or R2, or a gap if the code is not valid.
 */
static char DecodeBase(const char code, const std::string_view decoding_table) {
  if (IsMismatchToR1R2(code)) {
    return decoding_table.at(ToUnsigned(code - 'A'));
  }
  if (IsGapToR1R2(code)) {
    return decoding_table.at(ToUnsigned(code - 'a'));
  }
  return kDiscordantGapBase;
}

/**
 * @brief Decode a character code to a base in R1.
 * @param code Discordant character code
 * @return Decoded base in R1, or a gap if the code is not valid.
 */
static char DecodeR1(const char code) {
  return DecodeBase(code, kDecodingR1);
}

/**
 * @brief Decode a character code to a base in R2.
 * @param code Discordant character code
 * @return Decoded base in R2, or a gap if the code is not valid.
 */
static char DecodeR2(const char code) {
  return DecodeBase(code, kDecodingR2);
}

/**
 * @brief Helper function to deserialize a segment of the YC tag from an input string.
 * @param input_str String containing a segment of the YC tag
 * @return Vector of YcTagDuplexOp for the segment.
 */
static vec<YcTagDuplexOp> DeserializeSegment(const std::string& input_str, const std::string& read_name) {
  vec<YcTagDuplexOp> ops;
  u64 num_op_len = 0;
  static constexpr u64 kBase10 = 10;
  for (const char c : input_str) {
    if (c >= '0' && c <= '9') {
      // `c` is a digit
      if (num_op_len == 0) {
        num_op_len = static_cast<u64>(c - '0');
      } else {
        // accumulate the number
        num_op_len = num_op_len * kBase10 + static_cast<u64>(c - '0');
      }
    } else {
      if (num_op_len > 0) {
        // store the previous numeric op
        ops.emplace_back(true, 0, num_op_len);
        num_op_len = 0;
      }
      if (IsMismatchToR1R2(c)) {
        ops.emplace_back(false, c, 1);
      } else if (IsGapToR1R2(c)) {
        ops.emplace_back(false, c, 0);
      } else {
        throw error::Error(
            "Failed to deserialize YC tag for read named '{}': Invalid code '{}' in YC tag '{}'; YC tag might be "
            "malformed",
            read_name,
            c,
            input_str);
      }
    }
  }
  if (num_op_len > 0) {
    // store the final numeric op
    ops.emplace_back(true, 0, num_op_len);
  }
  return ops;
}

YcTag DeserializeYcTag(const std::string& input_str, const std::string& read_name) {
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

  // tag delimiter is not found, found once, or found >= 3x
  if (delimiter_pos1 == std::string::npos || delimiter_pos1 == delimiter_pos2 ||
      input_str.find_first_of(kTagDelimiters, delimiter_pos1 + 1) != delimiter_pos2) {
    throw error::Error(
        "Failed to deserialize YC tag for read named '{}': YC tag '{}' is not in the correct format; expected 3 "
        "components delimited by `+` or `-`",
        read_name,
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
      throw error::Error(
          "Failed to deserialize YC tag for read named '{}': YC tag '{}' is not in the correct format; left overhang "
          "length '{}' is not a number",
          read_name,
          input_str,
          left_overhang);
    }
    left_overhang_len = std::stoul(left_overhang);
  }
  if (!right_overhang.empty()) {
    if (right_overhang.find_first_not_of(kNumbers) != std::string::npos) {
      throw error::Error(
          "Failed to deserialize YC tag for read named '{}': YC tag '{}' is not in the correct format; right overhang "
          "length '{}' is not a number",
          read_name,
          input_str,
          right_overhang);
    }
    right_overhang_len = std::stoul(right_overhang);
  }

  return YcTag{left_overhang_len,
               DeserializeSegment(duplex_segment, read_name),
               right_overhang_len,
               input_str[delimiter_pos1] == kTagDelimiterR1,
               input_str[delimiter_pos2] == kTagDelimiterR1};
}

std::string GetYcTagString(const bam1_t* const record) {
  // extract YC tag string from BAM record aux fields
  const u8* const yc_tag_data = bam_aux_get(record, "YC");
  if (yc_tag_data == nullptr) {
    return "";  // return empty string if YC tag is not present
  }
  return std::string{bam_aux2Z(yc_tag_data)};
}

YcTag DeserializeYcTag(const bam1_t* const record) {
  const auto yc_tag_str = GetYcTagString(record);
  if (yc_tag_str.empty()) {
    return YcTag{};  // return empty YC tag
  }
  // extract hard clip lengths from CIGAR string
  u64 left_hard_clip = 0;
  u64 right_hard_clip = 0;
  if (const auto cigar_ops = record->core.n_cigar; cigar_ops > 1) {
    const u32* const cigar = bam_get_cigar(record);
    if (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP) {
      left_hard_clip = bam_cigar_oplen(cigar[0]);
    }
    if (bam_cigar_op(cigar[cigar_ops - 1]) == BAM_CHARD_CLIP) {
      right_hard_clip = bam_cigar_oplen(cigar[cigar_ops - 1]);
    }
  }
  // deserialize YC tag string and process it as necessary
  const std::string read_name = bam_get_qname(record);
  auto yc_tag = DeserializeYcTag(yc_tag_str, read_name);
  if (bam_is_rev(record)) {
    // if reverse-complement is required, it must be performed before trimming
    yc_tag.ReverseComplement(read_name);
  }
  if ((left_hard_clip > 0 || right_hard_clip > 0) && std::cmp_less(record->core.l_qseq, yc_tag.GetSequenceLength())) {
    // BAM record has hard clip(s), and the YC tag is longer than the consensus sequence.
    // Trim the YC tag based on the hard clips.
    yc_tag.Trim(left_hard_clip, right_hard_clip);
  }
  return yc_tag;
}

/**
 * @brief Helper function to reverse-complement a segment of the YC tag in place.
 * @param ops Vector of YcTagOp to be reverse-complemented
 */
static void ReverseComplementSegment(vec<YcTagDuplexOp>& ops, const std::string& read_name) {
  if (!ops.empty()) {
    // reverse the ops and complement any valid codes
    std::ranges::reverse(ops);
    for (auto& op : ops) {
      if (op.is_numeric) {
        continue;
      }
      if (IsMismatchToR1R2(op.code)) {
        op.code = kComplement[op.code - 'A'];
      } else if (IsGapToR1R2(op.code)) {
        op.code = kGapComplement[op.code - 'a'];
      } else {
        throw error::Error(
            "Failed to reverse complement YC tag for read named '{}': Invalid code '{}' in YC tag found during reverse "
            "complement operation",
            read_name,
            op.code);
      }
    }
  }
}

void YcTag::ReverseComplement(const std::string& read_name) {
  // swap left and right overhangs and their read assignments, but do not flip the read assignments (e.g. R1 -> R2)
  // because reverse-complemented R1 and R2 are still R1 and R2
  std::swap(left_overhang, right_overhang);
  std::swap(left_overhang_is_r1, right_overhang_is_r1);
  ReverseComplementSegment(duplex_ops, read_name);
}

u64 YcTag::GetSequenceLength() const {
  u64 length = left_overhang + right_overhang;
  for (const auto& op : duplex_ops) {
    length += op.length;
  }
  return length;
}

/**
 * @brief Helper function to trim a segment of the YC tag in place from the left.
 * @param ops YC tag ops of the segment to be trimmed
 * @param clip_len Length to trim
 * @return Remaining clip length after trimming; positive if clip length is not fully consumed, zero if fully consumed.
 */
static u64 LeftTrimSegment(vec<YcTagDuplexOp>& ops, u64 clip_len) {
  u64 i = 0;
  const auto ops_size = ops.size();
  while (i < ops_size && clip_len > 0) {
    if (clip_len < ops[i].length) {
      // trim the op by clip length
      ops[i].length -= clip_len;
      clip_len = 0;
      // break the loop because the op is not fully consumed
      break;
    }
    // consume the entire op
    clip_len -= ops[i].length;
    ops[i].length = 0;
    ++i;
  }
  // remove consumed leading ops
  ops.erase(ops.begin(), ops.begin() + static_cast<s32>(i));
  return clip_len;
}

/**
 * @brief Helper function to trim a segment of the YC tag in place from the right.
 * @param ops YC tag ops of the segment to be trimmed
 * @param clip_len Length to trim
 * @return Remaining clip length after trimming; positive if clip length is not fully consumed, zero if fully consumed.
 */
static u64 RightTrimSegment(vec<YcTagDuplexOp>& ops, u64 clip_len) {
  u64 i = ops.size();
  while (i > 0 && clip_len > 0) {
    const u64 idx = i - 1;
    if (clip_len < ops[idx].length) {
      // trim the op by clip length
      ops[idx].length -= clip_len;
      clip_len = 0;
      // break the loop because the op is not fully consumed
      break;
    }
    // consume the entire op
    clip_len -= ops[idx].length;
    ops[idx].length = 0;
    --i;
  }
  // remove consumed trailing ops
  ops.erase(ops.begin() + static_cast<s64>(i), ops.end());
  return clip_len;
}

void YcTag::Trim(u64 left_clip, u64 right_clip) {
  if (left_clip == 0 && right_clip == 0) {
    // no trimming needed
    return;
  }

  // consume left-clip from left to right
  if (left_clip < left_overhang) {
    // trim overhang and consume entire clip length
    left_overhang -= left_clip;
    left_clip = 0;
  } else {
    // consume entire overhang and reduce clip length
    left_clip -= left_overhang;
    left_overhang = 0;
  }
  if (left_clip > 0 && !duplex_ops.empty()) {
    left_clip = LeftTrimSegment(duplex_ops, left_clip);
  }
  if (left_clip < right_overhang) {
    // consume entire clip length
    right_overhang -= left_clip;
  } else {
    // consume entire overhang
    right_overhang = 0;
  }

  // consume right-clip from right to left
  if (right_clip < right_overhang) {
    // trim overhang and consume entire clip length
    right_overhang -= right_clip;
    right_clip = 0;
  } else {
    // consume entire overhang and reduce clip length
    right_clip -= right_overhang;
    right_overhang = 0;
  }
  if (right_clip > 0 && !duplex_ops.empty()) {
    right_clip = RightTrimSegment(duplex_ops, right_clip);
  }
  if (right_clip < left_overhang) {
    // trim overhang
    left_overhang -= right_clip;
  } else {
    // consume entire overhang
    left_overhang = 0;
  }
}

bool YcTag::IsDuplex() const {
  return !duplex_ops.empty();
}

bool YcTag::IsSimplex() const {
  return duplex_ops.empty();
}

u64 YcTag::CountDiscordantBase() const {
  u64 count{0};
  for (const auto& op : duplex_ops) {
    if (!op.is_numeric && DecodeR1(op.code) != DecodeR2(op.code)) {
      ++count;
    }
  }
  return count;
}

u64 YcTag::CountDuplexBase() const {
  return std::accumulate(duplex_ops.begin(), duplex_ops.end(), 0UL, [](const u64 sum, const YcTagDuplexOp& op) {
    return sum + op.length;
  });
}

vec<BaseType> YcTag::GetBaseTypes() const {
  using enum BaseType;
  vec<BaseType> base_types{};
  base_types.reserve(GetSequenceLength());
  if (left_overhang > 0) {
    base_types.insert(base_types.begin(), left_overhang, kSimplex);
  }
  for (const auto& op : duplex_ops) {
    if (op.length > 0) {
      if (op.is_numeric) {
        base_types.insert(base_types.end(), op.length, kConcordant);
      } else {
        // Each code can indicate either a concordant or discordant consensus base.
        // A consensus base is concordant if the decoded bases are the same between R1 and R2.
        // A consensus base is discordant if the decoded bases are different between R1 and R2.
        // Example codes: 'TtT'
        // There are 2 concordant 'T' bases in R1 and R2, but the consensus has non-T bases at these positions.
        // The gap code 't' is not assigned a base type because it has length 0.
        const auto base_type = DecodeR1(op.code) == DecodeR2(op.code) ? kConcordant : kDiscordant;
        base_types.insert(base_types.end(), op.length, base_type);
      }
    }
  }
  if (right_overhang > 0) {
    base_types.insert(base_types.end(), right_overhang, kSimplex);
  }
  return base_types;
}

std::string YcTag::ToString() const {
  std::string yc_tag_str{};
  if (left_overhang > 0) {
    // serialize left overhang only if it is not empty
    yc_tag_str += std::to_string(left_overhang);
  }
  yc_tag_str += left_overhang_is_r1 ? kTagDelimiterR1 : kTagDelimiterR2;
  for (const auto& [is_numeric, code, length] : duplex_ops) {
    if (is_numeric) {
      if (length > 0) {
        // serialize numeric op only if length > 0
        yc_tag_str += std::to_string(length);
      }
    } else {
      yc_tag_str += code;
    }
  }
  yc_tag_str += right_overhang_is_r1 ? kTagDelimiterR1 : kTagDelimiterR2;
  if (right_overhang > 0) {
    // serialize right overhang only if it is not empty
    yc_tag_str += std::to_string(right_overhang);
  }
  return yc_tag_str;
}

std::string EncodeBaseTypesToYcTag(const vec<BaseType>& base_types, const bool reverse, const std::string& read_name) {
  using enum BaseType;
  u64 left_overhang = 0;
  std::string duplex_segment{};
  u64 num_simplex = 0;
  u64 num_concordant = 0;

  // lamda function to process a given base type
  auto process_type =
      [&left_overhang, &duplex_segment, &num_simplex, &num_concordant, &read_name](const BaseType& base_type) {
        if (base_type == kDiscordant) {
          if (num_simplex > 0) {
            left_overhang = std::exchange(num_simplex, 0);
          }
          if (num_concordant > 0) {
            duplex_segment.append(std::to_string(std::exchange(num_concordant, 0)));
          }
          duplex_segment += kDiscordantChar;
        } else if (base_type == kConcordant) {
          if (num_simplex > 0) {
            left_overhang = std::exchange(num_simplex, 0);
          }
          ++num_concordant;
        } else if (base_type == kSimplex) {
          ++num_simplex;
        } else {
          throw BaseTypeError(read_name, base_type);
        }
      };

  // process base types in the specified order
  if (reverse) {
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
    left_overhang = std::exchange(num_simplex, 0);
  }

  const std::string left_overhang_str = left_overhang > 0 ? std::to_string(left_overhang) : "";
  const std::string right_overhang_str = num_simplex > 0 ? std::to_string(num_simplex) : "";
  // Use '+' for both overhang because the base types are identical regardless of the overhang assignments.
  return fmt::format("{}+{}+{}", left_overhang_str, duplex_segment, right_overhang_str);
}

std::string EncodeBaseTypesToYcTag(const vec<BaseType>& base_types, const std::string& read_name) {
  return EncodeBaseTypesToYcTag(base_types, false, read_name);
}

DecodedReads Decode(const std::string& consensus_seq,
                    const std::string& consensus_qual,
                    const YcTag& yc_tag,
                    const std::string& read_name) {
  const auto seq_len = consensus_seq.length();
  if (seq_len != consensus_qual.length()) {
    throw error::Error(
        "Failed to decode read named '{}': Consensus sequence indicates length of '{}' but consensus qualities has "
        "length of '{}'; this indicates a malformed consensus read from an upstream process",
        read_name,
        seq_len,
        consensus_qual.length());
  }
  const auto yc_tag_len = yc_tag.GetSequenceLength();
  if (seq_len != yc_tag_len) {
    throw error::Error(
        "Failed to decode read named '{}': Consensus sequence indicates length of '{}' but YC tag has length of '{}'; "
        "this indicates a malformed consensus read from an upstream process",
        read_name,
        seq_len,
        yc_tag_len);
  }
  // Both output reads would be in the same orientation as the consensus sequence.
  std::string r1_seq{};
  std::string r2_seq{};
  std::string r1_qual{};
  std::string r2_qual{};
  using enum BaseType;
  vec<BaseType> r1_base_types{};
  vec<BaseType> r2_base_types{};

  // lambda function to append simplex segment to either R1 or R2
  auto append_simplex =
      [&consensus_seq, &consensus_qual, &r1_seq, &r2_seq, &r1_qual, &r2_qual, &r1_base_types, &r2_base_types](
          const bool is_r1, const u64 start, const u64 len) {
        if (is_r1) {
          r1_seq += consensus_seq.substr(start, len);
          r1_qual += consensus_qual.substr(start, len);
          r1_base_types.insert(r1_base_types.end(), len, kSimplex);
        } else {
          r2_seq += consensus_seq.substr(start, len);
          r2_qual += consensus_qual.substr(start, len);
          r2_base_types.insert(r2_base_types.end(), len, kSimplex);
        }
      };

  // process the left overhang
  if (yc_tag.left_overhang > 0) {
    append_simplex(yc_tag.left_overhang_is_r1, 0, yc_tag.left_overhang);
  }

  // process the duplex segment
  u64 idx = yc_tag.left_overhang;
  for (const auto& [is_numeric, code, length] : yc_tag.duplex_ops) {
    if (is_numeric) {
      const auto seq = consensus_seq.substr(idx, length);
      r1_seq += seq;
      r2_seq += seq;
      const auto qual = consensus_qual.substr(idx, length);
      r1_qual += qual;
      r2_qual += qual;
      r1_base_types.insert(r1_base_types.end(), length, kConcordant);
      r2_base_types.insert(r2_base_types.end(), length, kConcordant);
      idx += length;
      continue;
    }
    // Not numeric, this is a discordant base.

    // If consensus base is not a gap to R1/R2, extract the consensus base quality
    const auto qual = (length == 1 && idx < consensus_qual.length()) ? std::string(1, consensus_qual.at(idx))
                                                                     : kConsensusGapBaseQuality;
    const auto r1_base = DecodeR1(code);
    const auto r2_base = DecodeR2(code);
    const bool same_base = r1_base == r2_base;

    // lambda function to process a decoded base
    auto process_base = [&qual, &same_base](
                            const char base, std::string& seq, std::string& qual_str, vec<BaseType>& base_types) {
      if (base == kDiscordantGapBase) {
        if (!base_types.empty()) {
          // The gap does not add a base (or quality score or base type) to R1, but a base is added to R2 instead.
          // To indicate this for R1, modify the base type at the previous position in R1 to "discordant".
          // A distinction between discordant and concordant indels is crucial to accurate variant calling.
          base_types.back() = kDiscordant;
        }
      } else {
        seq += base;
        qual_str += qual;
        // If the decoded base is identical for R1 and R2 but not the same as the consensus base, then
        // the decoded base type is still considered as "concordant" for both R1 and R2.
        base_types.push_back(same_base ? kConcordant : kDiscordant);
      }
    };

    process_base(r1_base, r1_seq, r1_qual, r1_base_types);
    process_base(r2_base, r2_seq, r2_qual, r2_base_types);
    idx += length;
  }

  // process the right overhang
  if (idx < seq_len) {
    append_simplex(yc_tag.right_overhang_is_r1, idx, yc_tag.right_overhang);
  }
  return DecodedReads{r1_seq, r2_seq, r1_qual, r2_qual, r1_base_types, r2_base_types};
}

DecodedReads Decode(const std::string& consensus_seq,
                    const std::string& consensus_qual,
                    const std::string& yc_tag_str,
                    const bool rev_comp,
                    const u64 left_clip,
                    const u64 right_clip,
                    const std::string& read_name) {
  YcTag yc_tag = DeserializeYcTag(yc_tag_str, read_name);
  if (rev_comp) {
    // If reverse-complement is required, it must be done before trimming
    yc_tag.ReverseComplement(read_name);
  }
  if ((left_clip > 0 || right_clip > 0) && consensus_seq.size() < yc_tag.GetSequenceLength()) {
    // There are hard clips, and the YC tag is longer than the consensus sequence.
    // Trim the YC tag based on the hard clips.
    yc_tag.Trim(left_clip, right_clip);
  }
  return Decode(consensus_seq, consensus_qual, yc_tag, read_name);
}

DecodedReads Decode(const std::string& consensus_seq,
                    const std::string& consensus_qual,
                    const std::string& yc_tag_str,
                    const bool rev_comp,
                    const std::string& read_name) {
  return Decode(consensus_seq, consensus_qual, yc_tag_str, rev_comp, 0, 0, read_name);
}

DecodedReads Decode(const std::string& consensus_seq,
                    const std::string& consensus_qual,
                    const std::string& yc_tag_str,
                    const std::string& read_name) {
  return Decode(consensus_seq, consensus_qual, yc_tag_str, false, 0, 0, read_name);
}

/**
 * @brief Helper function to add a CIGAR operation to a CIGAR vector, merging with the last operation if they are the
 * same.
 * @param cigar Vector of CIGAR operations
 * @param op CIGAR operation
 * @param len Length of the CIGAR operation
 */
static void AddCigarOp(vec<HtslibCigarOp>& cigar, u32 op, u64 len) {
  if (len > 0) {
    if (!cigar.empty() && cigar.back().op == op) {
      cigar.back().len += len;
    } else {
      cigar.emplace_back(op, len);
    }
  }
}

/**
 * @brief Helper function to handle a gap in R1 or R2 when updating the CIGAR operations.
 * @param cigar CIGAR vector to update
 * @param cigar_op CIGAR operation code
 * @param cigar_op_consumes_ref Boolean indicating if the CIGAR operation consumes reference bases
 * @param yc_op_len Length of the YC operation
 * @param decoded_base_is_gap Boolean indicating if the current decoded base is a gap
 * @param adjust_ref_pos Boolean reference to indicate if reference position adjustment is needed
 * @param ref_pos Reference position to update
 */
static void HandleGap(vec<HtslibCigarOp>& cigar,
                      const u32 cigar_op,
                      const bool cigar_op_consumes_ref,
                      const u64 yc_op_len,
                      const bool decoded_base_is_gap,
                      bool& adjust_ref_pos,
                      u64& ref_pos) {
  if (decoded_base_is_gap) {
    if (cigar_op_consumes_ref) {
      if (adjust_ref_pos) {
        ++ref_pos;
      }
      AddCigarOp(cigar, BAM_CDEL, yc_op_len);
    }
  } else {
    adjust_ref_pos = false;
    AddCigarOp(cigar, cigar_op, yc_op_len);
  }
}

/**
 * @brief Helper function to update the alignment information for the duplex segment of the YC tag.
 * @param yc_tag YC tag object
 * @param yc_op_read_start Starting position of the read in the YC tag
 * @param cigar_infos Vector of CIGAR align info
 * @param cigar_infos_idx Current index in the CIGAR align info
 * @param r1_cigar Vector to store updated CIGAR operations for R1
 * @param r2_cigar Vector to store updated CIGAR operations for R2
 * @param r1_ref_pos Reference position for R1
 * @param r2_ref_pos Reference position for R2
 */
static void UpdateDuplexSegmentAlignmentInfo(const YcTag& yc_tag,
                                             u64& yc_op_read_start,
                                             const vec<HtslibCigarAlignInfo>& cigar_infos,
                                             u64& cigar_infos_idx,
                                             vec<HtslibCigarOp>& r1_cigar,
                                             vec<HtslibCigarOp>& r2_cigar,
                                             u64& r1_ref_pos,
                                             u64& r2_ref_pos) {
  // If left overhang is empty, then the duplex segment starts at the beginning of the read.
  // Adjust reference positions if the leading YC tag ops decode to gap(s) in R1/R2 with respect to the consensus.
  bool adjust_r1_ref_pos = !yc_tag.left_overhang_is_r1 || yc_tag.left_overhang == 0;
  bool adjust_r2_ref_pos = yc_tag.left_overhang_is_r1 || yc_tag.left_overhang == 0;

  for (const auto& [yc_op_is_numeric, yc_op_code, yc_op_len] : yc_tag.duplex_ops) {
    const char decoded_r1_base = DecodeR1(yc_op_code);
    const char decoded_r2_base = DecodeR2(yc_op_code);

    if (yc_op_len == 0) {
      // consensus base is a gap with respect to R1/R2
      if (decoded_r1_base != kDiscordantGapBase) {
        adjust_r1_ref_pos = false;
        AddCigarOp(r1_cigar, BAM_CINS, 1);
      }
      if (decoded_r2_base != kDiscordantGapBase) {
        adjust_r2_ref_pos = false;
        AddCigarOp(r2_cigar, BAM_CINS, 1);
      }
      continue;
    }

    // consensus base is not a gap with respect to R1/R2
    const auto yc_op_read_end = yc_op_read_start + yc_op_len;
    for (; cigar_infos_idx < cigar_infos.size(); ++cigar_infos_idx) {
      const auto& [cigar_op, read_pos, ref_pos, read_span, ref_span] = cigar_infos.at(cigar_infos_idx);
      if (read_pos >= yc_op_read_end) {
        // CIGAR op starts after the duplex segment
        break;
      }

      const bool cigar_op_consumes_read = ConsumesRead(cigar_op);
      const bool cigar_op_consumes_ref = ConsumesReference(cigar_op);

      if (cigar_op_consumes_ref && !cigar_op_consumes_read) {
        // CIGAR op consumes reference base(s) but not read base(s)
        adjust_r1_ref_pos = false;
        adjust_r2_ref_pos = false;
        AddCigarOp(r1_cigar, cigar_op, ref_span);
        AddCigarOp(r2_cigar, cigar_op, ref_span);
      } else if (yc_op_is_numeric) {
        // CIGAR op consumes read base(s)
        adjust_r1_ref_pos = false;
        adjust_r2_ref_pos = false;
        const u64 span = std::min(read_pos + read_span, yc_op_read_end) - std::max(read_pos, yc_op_read_start);
        AddCigarOp(r1_cigar, cigar_op, span);
        AddCigarOp(r2_cigar, cigar_op, span);
      } else {
        // YC tag op is not numeric, i.e., a discordant base
        HandleGap(r1_cigar,
                  cigar_op,
                  cigar_op_consumes_ref,
                  yc_op_len,
                  decoded_r1_base == kDiscordantGapBase,
                  adjust_r1_ref_pos,
                  r1_ref_pos);
        HandleGap(r2_cigar,
                  cigar_op,
                  cigar_op_consumes_ref,
                  yc_op_len,
                  decoded_r2_base == kDiscordantGapBase,
                  adjust_r2_ref_pos,
                  r2_ref_pos);
      }

      if (read_pos + read_span > yc_op_read_end) {
        // CIGAR op spans beyond the duplex segment, do not look at the next CIGAR op
        break;
      }
    }
    yc_op_read_start = yc_op_read_end;
  }
}

/**
 * @brief Helper function to reconcile neighboring deletion and insertion CIGAR ops.
 * It removes trailing deletions, converts leading insertions to soft-clips,
 * and reconciles neighboring deletion and insertion ops into match ops when appropriate.
 * @param ops Vector of CIGAR operations to be processed
 * @return Processed vector of CIGAR operations
 */
static vec<HtslibCigarOp> ProcessIndelCigarOps(const vec<HtslibCigarOp>& ops) {
  if (ops.empty()) {
    return {};
  }

  vec<HtslibCigarOp> new_ops;
  auto it = ops.begin();

  // Lambda function to add a CIGAR ops to the `new_ops` vector, merging with the last op if they are the same.
  auto add_op = [&new_ops](const u32 op, const u32 len) {
    if (len == 0) {
      return;
    }
    if (!new_ops.empty() && new_ops.back().op == op) {
      new_ops.back().len += len;
    } else {
      new_ops.emplace_back(op, len);
    }
  };

  while (it != ops.end()) {
    if (it->op == BAM_CDEL || it->op == BAM_CINS) {
      // Tally consecutive deletion and insertion op lengths
      u32 del_len = 0;
      u32 ins_len = 0;
      for (; it != ops.end() && (it->op == BAM_CDEL || it->op == BAM_CINS); ++it) {
        if (it->op == BAM_CDEL) {
          del_len += it->len;
        } else {
          ins_len += it->len;
        }
      }
      // Neighboring indel ops are converted to a single match op.
      // Any remaining indel lengths are added as a separate op after the match op.
      // This is okay because the match op represents alignment match, which can include sequence mismatches.
      const u32 match_len = std::min(del_len, ins_len);
      u32 ins_remain = ins_len > del_len ? ins_len - del_len : 0;
      const u32 del_remain = del_len > ins_len ? del_len - ins_len : 0;

      if (ins_remain > 0 && new_ops.empty()) {
        // leading insertion op is converted to soft-clip op
        add_op(BAM_CSOFT_CLIP, ins_remain);
        ins_remain = 0;
      }
      add_op(BAM_CMATCH, match_len);
      add_op(BAM_CINS, ins_remain);
      if (del_remain > 0 && !new_ops.empty()) {
        // leading deletion ops are removed
        add_op(BAM_CDEL, del_remain);
      }
      continue;
    }
    add_op(it->op, it->len);
    ++it;
  }

  // Remove trailing deletion op and convert trailing insertion op to soft-clip
  if (!new_ops.empty()) {
    const auto& last_op = new_ops.back().op;
    if (last_op == BAM_CDEL) {
      // example: 10M2D -> 10M
      new_ops.pop_back();
    } else if (last_op == BAM_CINS) {
      // example: 10M2I -> 10M2S
      new_ops.back().op = BAM_CSOFT_CLIP;
    }
  }
  return new_ops;
}

/**
 * @brief Helper function to extract alignment information for decoded alignments.
 * @param yc_tag YC tag
 * @param record BAM record
 * @param decoded_alns Decoded alignments to be updated with alignment information
 */
static void UpdateAlignmentInfo(const YcTag& yc_tag, const bam1_t* const record, DecodedAlignments& decoded_alns) {
  const auto cigar_infos = GetCigarAlignInfos(record);
  if (cigar_infos.empty()) {
    return;
  }

  u64 cigar_info_idx = 0;
  u64 yc_tag_op_read_start = 0;

  const bool has_duplex = !yc_tag.duplex_ops.empty();
  u64 r1_ref_pos = ToUnsigned(record->core.pos);
  u64 r2_ref_pos = r1_ref_pos;
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

  if (has_duplex) {
    // process CIGAR op(s) within the duplex segment
    UpdateDuplexSegmentAlignmentInfo(
        yc_tag, yc_tag_op_read_start, cigar_infos, cigar_info_idx, r1_cigar, r2_cigar, r1_ref_pos, r2_ref_pos);
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

  decoded_alns.r1_ref_pos = decoded_alns.r1_seq.empty() ? 0 : r1_ref_pos;
  decoded_alns.r2_ref_pos = decoded_alns.r2_seq.empty() ? 0 : r2_ref_pos;
  decoded_alns.r1_cigar = ProcessIndelCigarOps(r1_cigar);
  decoded_alns.r2_cigar = ProcessIndelCigarOps(r2_cigar);
}

DecodedAlignments Decode(const bam1_t* const record) {
  const auto yc_tag = DeserializeYcTag(record);
  if (yc_tag.GetSequenceLength() == 0) {
    return DecodedAlignments{};
  }
  // extract consensus sequence and base quality
  std::string consensus_seq;
  std::string consensus_qual;
  const auto seq_len = record->core.l_qseq;
  if (seq_len <= 0) {
    throw error::Error(
        "Failed to decode read named '{}': Invalid sequence length of '{}'; this indicates a malformed consensus read "
        "from an upstream process",
        bam_get_qname(record),
        seq_len);
  }
  const auto seq_len_u = static_cast<size_t>(seq_len);
  consensus_seq.resize(seq_len_u);
  consensus_qual.resize(seq_len_u);
  const auto* const seq_data = bam_get_seq(record);
  const auto* const qual_data = bam_get_qual(record);
  for (size_t i = 0; i < seq_len_u; ++i) {
    consensus_seq[i] = seq_nt16_str[bam_seqi(seq_data, i)];
    consensus_qual[i] = static_cast<char>(qual_data[i] + kPhred33Offset);
  }
  // decode the consensus sequence using the YC tag
  const auto [r1_seq, r2_seq, r1_qual, r2_qual, r1_base_types, r2_base_types] =
      Decode(consensus_seq, consensus_qual, yc_tag, bam_get_qname(record));
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

/**
 * @brief Helper function to set up a new BAM record based on the source BAM record.
 * @param src_record Source BAM record
 * @param new_record New BAM record to be filled
 * @param seq Decoded sequence
 * @param qual Decoded quality
 * @param ref_pos Reference position for the decoded read
 * @param qname Read name
 * @param cigar_ops CIGAR operations for the decoded read
 * @param yc_tag YC tag string for the decoded read
 */
static void DecodeToBamRecordsHelper(const bam1_t* const src_record,
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
  const u8* const aux_first = bam_aux_first(src_record);
  const u8* aux = aux_first;
  while (aux != nullptr) {
    const u8* const next_aux = bam_aux_next(src_record, aux);
    // Determine the length of the auxiliary data
    const auto length = static_cast<s32>((next_aux != nullptr) ? next_aux - aux : aux_first + aux_data_len - aux);
    // Subtract 3 from the aux data length to account for the tag (2 bytes) and the type (1 byte)
    static constexpr s32 kAuxLengthOffset = 3;
    // Append the auxiliary data to the destination record
    bam_aux_append(new_record.get(), bam_aux_tag(aux), bam_aux_type(aux), length - kAuxLengthOffset, aux + 1);
    aux = next_aux;
  }
  bam_aux_update_str(new_record.get(), "YC", static_cast<s32>(yc_tag.length() + 1), yc_tag.c_str());
}

std::variant<std::monostate, io::Bam1Ptr, std::pair<io::Bam1Ptr, io::Bam1Ptr>> DecodeToBamRecords(
    const bam1_t* const record) {
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
                             EncodeBaseTypesToYcTag(decoded_alns.r1_base_types, rev_comp, qname));
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
                             EncodeBaseTypesToYcTag(decoded_alns.r2_base_types, rev_comp, qname));
    if (!has_r1) {
      return std::move(r2_record);
    }
  }
  return std::pair{std::move(r1_record), std::move(r2_record)};
}

/**
 * Ensures that the length of the read sequence matches the size of the provided base types vector.
 *
 * @param read Pointer to a BAM record (bam1_t) representing the read.
 * @param base_types A vector of yc_decode::BaseType representing the base types of the read sequence.
 *
 * @throws error::Error If the size of the base_types vector does not match the length of the read sequence.
 */
void RequireConsistentReadLength(const bam1_t* const read, const vec<yc_decode::BaseType>& base_types) {
  if (static_cast<size_t>(read->core.l_qseq) != base_types.size()) {
    throw error::Error(
        "Failed to process base types in read named '{}': Read sequence indicates length of '{}' but base types has "
        "length of '{}'; YC tag might be malformed",
        bam_get_qname(read),
        static_cast<size_t>(read->core.l_qseq),
        base_types.size());
  }
}

std::runtime_error BaseTypeError(const std::string& read_name, const BaseType& invalid_base_type) {
  return error::Error(
      "Failed to process base type in read named '{}': Invalid base type enum value '{}'; YC tag might be malformed or "
      "there is an issue encoding base types in the YC tag format",
      read_name,
      enum_util::FormatEnumName(invalid_base_type));
}

}  // namespace xoos::yc_decode
