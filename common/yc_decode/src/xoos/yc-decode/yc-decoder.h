#pragma once

#include <string>
#include <variant>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "cigar-utils.h"

namespace xoos::yc_decode {

// serves as a replacement for legacy base quality values 0,18,93 (or 5,22,39)
enum class BaseType {
  kDiscordant,
  kSimplex,
  kConcordant
};

// represents an operator within the duplex segment of the YC tag
struct YcTagDuplexOp {
  bool is_numeric{false};  // flag whether this op is numeric
  char code{0};            // Encoded character for discordant base; not required for numeric op
  u64 length{0};           // op length; numeric ops have variable lengths, all discordant ops have length of 1
  auto operator<=>(const YcTagDuplexOp&) const = default;
};

// information extracted from a YC tag string
struct YcTag {
  // member variables
  u64 left_overhang{0};              // length of left overhang
  vec<YcTagDuplexOp> duplex_ops{};   // duplex segment ops
  u64 right_overhang{0};             // length of right overhang
  bool left_overhang_is_r1{false};   // true if left overhang is R1
  bool right_overhang_is_r1{false};  // true if right overhang is R1
  // member functions
  u64 GetSequenceLength() const;
  void ReverseComplement();
  void Trim(u64 left_clip, u64 right_clip);
  bool IsDuplex() const;
  bool IsSimplex() const;
  u64 CountDiscordantBase() const;
  u64 CountDuplexBase() const;
  vec<BaseType> GetBaseTypes() const;
  std::string ToString() const;
  auto operator<=>(const YcTag&) const = default;
};

struct DecodedReads {
  std::string r1_seq{};
  std::string r2_seq{};
  std::string r1_qual{};
  std::string r2_qual{};
  vec<BaseType> r1_base_types{};
  vec<BaseType> r2_base_types{};
  auto operator<=>(const DecodedReads&) const = default;
};

struct DecodedAlignments : DecodedReads {
  u64 r1_ref_pos{0};              // 0-based start position in the reference for R1
  u64 r2_ref_pos{0};              // 0-based start position in the reference for R2
  vec<HtslibCigarOp> r1_cigar{};  // htslib CIGAR ops for R1
  vec<HtslibCigarOp> r2_cigar{};  // htslib CIGAR ops for R2
  auto operator<=>(const DecodedAlignments&) const = default;
};

constexpr u8 kPhred33Offset = 33;

YcTag DeserializeYcTag(const std::string& input_str);
YcTag DeserializeYcTag(const bam1_t* record);
std::string EncodeBaseTypesToYcTag(const vec<BaseType>& base_types, bool reverse = false);
DecodedReads Decode(const std::string& consensus_seq, const std::string& consensus_qual, const YcTag& yc_tag);
DecodedReads Decode(const std::string& consensus_seq,
                    const std::string& consensus_qual,
                    const std::string& yc_tag_str,
                    bool rev_comp = false,
                    u64 left_clip = 0,
                    u64 right_clip = 0);
DecodedAlignments Decode(const bam1_t* record);
std::variant<std::monostate, io::Bam1Ptr, std::pair<io::Bam1Ptr, io::Bam1Ptr>> DecodeToBamRecords(const bam1_t* record);

std::runtime_error BaseTypeError(const std::string& read_name);

}  // namespace xoos::yc_decode
