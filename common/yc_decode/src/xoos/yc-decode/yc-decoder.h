#pragma once

#include <string>
#include <variant>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "cigar-utils.h"

namespace xoos::yc_decode {

/**
 * Synopsis:
 * A duplex consensus sequence is derived from its two constituent reads (R1, R2). The YC tag is used to encode the
 * differences between the consensus sequence and R1/R2. This library provides functions to decode YC tags and
 * reconstruct R1 and R2 sequences, qualities, and alignments.
 */

/**
 * @brief Base types for duplex consensus read.
 * @details This serves as a replacement for legacy base quality values 0,18,93 (or 5,22,39).
 */
enum class BaseType {
  kDiscordant,
  kSimplex,
  kConcordant
};

/**
 * @brief An operator within the duplex segment of a YC tag.
 * @details An operator can be either a number or character code. Numeric ops have variable lengths and their codes are
 * always ignored. Non-numeric ops with UPPER case character code have a length of 1. Non-numeric ops with lower case
 * character code have a length of 0.
 */
struct YcTagDuplexOp {
  bool is_numeric{false};
  // encoded character for non-numeric op
  char code{0};
  // op length in the consensus sequence
  u64 length{0};
  bool operator==(const YcTagDuplexOp&) const = default;
};

/**
 * @brief Information extracted from a YC tag string.
 */
struct YcTag {
  // length of left overhang
  u64 left_overhang{0};
  // duplex segment ops
  vec<YcTagDuplexOp> duplex_ops{};
  // length of right overhang
  u64 right_overhang{0};
  // flag whether left overhang is R1, false for R2
  bool left_overhang_is_r1{false};
  // flag whether right overhang is R1, false for R2
  bool right_overhang_is_r1{false};

  /**
   * @brief Return the total sequence length of the YC tag.
   * @return Total sequence length
   */
  u64 GetSequenceLength() const;

  /**
   * @brief Reverse-complement YC tag in place.
   */
  void ReverseComplement(const std::string& read_name);

  /**
   * @brief Trim both ends of YC tag in place.
   * @param left_clip Left clip length
   * @param right_clip Right clip length
   * @details If the clip length(s) were longer than the YC tag, the YC tag would be completely trimmed.
   */
  void Trim(u64 left_clip, u64 right_clip);

  /**
   * @brief Check whether the consensus read is duplex.
   * @return True if the read is duplex.
   */
  bool IsDuplex() const;

  /**
   * @brief Check whether the consensus read is simplex. Special case: an empty YC tag is considered simplex.
   * @return True if the read is simplex.
   */
  bool IsSimplex() const;

  /**
   * @brief Count the number of discordant bases between R1 and R2 in the duplex segment.
   * @return Number of discordant bases
   */
  u64 CountDiscordantBase() const;

  /**
   * @brief Count the number of bases in the duplex segment, including both concordant and discordant bases.
   * @return Number of duplex bases
   */
  u64 CountDuplexBase() const;

  /**
   * @brief Extract the base type at each position represented by the YC tag.
   * @return Vector of base type enums
   */
  vec<BaseType> GetBaseTypes() const;

  /**
   * @brief Serialize the YC tag to string format.
   * @return string format of the YC tag
   */
  std::string ToString() const;

  bool operator==(const YcTag&) const = default;
};

/**
 * @brief Decoded reads for a duplex consensus sequence based on its YC tag.
 */
struct DecodedReads {
  std::string r1_seq{};
  std::string r2_seq{};
  std::string r1_qual{};
  std::string r2_qual{};
  vec<BaseType> r1_base_types{};
  vec<BaseType> r2_base_types{};
  auto operator<=>(const DecodedReads&) const = default;
};

/**
 * @brief Decoded alignments for a duplex consensus sequence based on its YC tag.
 */
struct DecodedAlignments : DecodedReads {
  // 0-based start position in the reference for R1
  u64 r1_ref_pos{0};
  // 0-based start position in the reference for R2
  u64 r2_ref_pos{0};
  // htslib CIGAR ops for R1
  vec<HtslibCigarOp> r1_cigar{};
  // htslib CIGAR ops for R2
  vec<HtslibCigarOp> r2_cigar{};

  auto operator<=>(const DecodedAlignments&) const = default;
};

// Offset for converting base quality's integer value to PHRED33
constexpr u8 kPhred33Offset = 33;

/**
 * @brief Deserialize YC tag from input string.
 * @param input_str input string containing YC tag, e.g. "123-45M67-89"
 * @param read_name string containing read name for descriptive error throwing
 * @return YC tag object
 */
YcTag DeserializeYcTag(const std::string& input_str, const std::string& read_name);

/**
 * @brief Extract YC tag value string directly from a BAM record.
 * @param record BAM record pointer
 * @return YC tag value string
 */
std::string GetYcTagString(const bam1_t* record);

/**
 * @brief Deserialize YC tag from BAM record, reverse-complement and trim it as necessary.
 * @note If the sequence is aligned to the reverse strand, then the YC tag would be reverse-complemented.
 * @note YC tag that is longer than the sequence may be trimmed based on left and right hard clips.
 * @post If the YC tag is absent or empty, an empty YC tag object would be returned.
 * @param record BAM record pointer
 * @return YC tag object
 * @throws error::Error If the YC tag is present but malformed (e.g. invalid delimiter/overhang format or invalid duplex
 * codes).
 * @throws error::Error If reverse-complementing the YC tag encounters an unexpected code.
 */
YcTag DeserializeYcTag(const bam1_t* record);

/**
 * @brief Encode base types to YC tag string format, with option to reverse the order.
 * @param base_types Base types
 * @param reverse Flag whether to reverse the order of base types
 * @param read_name Read name for descriptive error throwing
 * @return YC tag string
 */
std::string EncodeBaseTypesToYcTag(const vec<BaseType>& base_types, bool reverse, const std::string& read_name);

/**
 * @brief Encode base types to YC tag string format.
 * @param base_types Base types
 * @param read_name Read name for descriptive error throwing
 * @return YC tag string
 */
std::string EncodeBaseTypesToYcTag(const vec<BaseType>& base_types, const std::string& read_name);

/**
 * @brief Decode the YC tag string for a consensus sequence into two constituent reads.
 * @pre Consensus sequence, base quality, and YC tag must be in the same orientation.
 * @pre Consensus sequence, base quality, and the underlying sequence represented by the YC tag must have equal length.
 * @post Decoded reads will be in the same orientation as the consensus sequence.
 * @post Left and right overhangs, if present, are assigned to R1 or R2 according to the YC tag delimiters (+/-).
 * @param consensus_seq Consensus sequence
 * @param consensus_qual Consensus quality
 * @param yc_tag YC tag struct
 * @param read_name Read name for descriptive error throwing
 * @return Decoded reads
 * @throws error::Error if consensus sequence and base quality have inconsistent length
 * @throws error::Error if consensus sequence and the underlying sequence represented by the YC tag have inconsistent
 * length.
 * @see Decode(const bam1_t* record) for decoding directly from BAM record. Hard clips and reverse-complement will be
 * handled automatically.
 */
DecodedReads Decode(const std::string& consensus_seq,
                    const std::string& consensus_qual,
                    const YcTag& yc_tag,
                    const std::string& read_name);

/**
 * @brief Decode the YC tag string for a consensus sequence into two constituent reads.
 * @pre Consensus sequence and base quality must be in the same orientation.
 * @pre Consensus sequence and base quality must have the same length.
 * @note YC tag in the opposite orientation will be reverse-complemented.
 * @note YC tag that is longer than the consensus sequence may be trimmed based on left/right clips.
 * @post Decoded reads will be in the same orientation as the consensus sequence.
 * @post Left and right overhangs, if present, are assigned to R1 or R2 according to the YC tag delimiters (+/-).
 * @param consensus_seq Consensus sequence
 * @param consensus_qual Consensus quality
 * @param yc_tag_str YC tag string
 * @param rev_comp Flag whether YC tag needs to be reverse-complemented
 * @param left_clip Left clip length for trimming YC tag
 * @param right_clip Right clip length for trimming YC tag
 * @param read_name Read name for descriptive error throwing
 * @return Decoded reads
 * @throws error::Error if consensus sequence and base quality have inconsistent length
 * @throws error::Error if consensus sequence and the underlying sequence represented by the YC tag still have
 * inconsistent length after trimming YC tag based on hard-clip length(s).
 * @see Decode(const bam1_t* record) for decoding directly from BAM record. Hard clips and reverse-complement will be
 * handled automatically.
 */
DecodedReads Decode(const std::string& consensus_seq,
                    const std::string& consensus_qual,
                    const std::string& yc_tag_str,
                    bool rev_comp,
                    u64 left_clip,
                    u64 right_clip,
                    const std::string& read_name);

/**
 * @brief Decode the YC tag string for a consensus sequence into two constituent reads.
 * @pre Consensus sequence and base quality must be in the same orientation.
 * @pre Consensus sequence, base quality, and the underlying sequence represented by the YC tag must have the same
 * length.
 * @note YC tag in the opposite orientation will be reverse-complemented.
 * @post Decoded reads will be in the same orientation as the consensus sequence.
 * @post Left and right overhangs, if present, are assigned to R1 or R2 according to the YC tag delimiters (+/-).*
 * @param consensus_seq Consensus sequence
 * @param consensus_qual Consensus quality
 * @param yc_tag_str YC tag string
 * @param rev_comp Flag whether YC tag needs to be reverse-complemented
 * @param read_name Read name for descriptive error throwing
 * @return Decoded reads
 * @throws error::Error if consensus sequence and base quality have inconsistent length
 * @throws error::Error if consensus sequence and the underlying sequence represented by the YC tag have inconsistent
 * length.
 * @throws error::Error if the YC tag string is malformed or cannot be parsed (e.g. invalid YC tag format or codes).
 * @throws error::Error if reverse-complementing the YC tag fails due to invalid or unsupported YC tag codes.
 * @see Decode(const bam1_t* record) for decoding directly from BAM record. Hard clips and reverse-complement will be
 * handled automatically.
 */
DecodedReads Decode(const std::string& consensus_seq,
                    const std::string& consensus_qual,
                    const std::string& yc_tag_str,
                    bool rev_comp,
                    const std::string& read_name);

/**
 * @brief Decode the YC tag string for a consensus sequence into two constituent reads.
 * @pre Consensus sequence, base quality, and YC tag must be in the same orientation.
 * @pre Consensus sequence, base quality, and the underlying sequence represented by the YC tag must have the same
 * length.
 * @post Decoded reads will be in the same orientation as the consensus sequence.
 * @post Left and right overhangs, if present, are assigned to R1 or R2 according to the YC tag delimiters (+/-).
 * @param consensus_seq Consensus sequence
 * @param consensus_qual Consensus quality
 * @param yc_tag_str YC tag string
 * @param read_name Read name for descriptive error throwing
 * @return Decoded reads
 * @throws error::Error if consensus sequence and base quality have inconsistent length
 * @throws error::Error if consensus sequence and the underlying sequence represented by the YC tag have inconsistent
 * length.
 * @throws error::Error if the YC tag string is malformed or cannot be parsed (e.g. invalid YC tag format or codes).
 * @see Decode(const bam1_t* record) for decoding directly from BAM record. Hard clips and reverse-complement will be
 * handled automatically.
 */
DecodedReads Decode(const std::string& consensus_seq,
                    const std::string& consensus_qual,
                    const std::string& yc_tag_str,
                    const std::string& read_name);

/**
 * @brief Decode a BAM record using its YC tag.
 * @pre Input BAM record must have a non-empty sequence.
 * @note If the sequence is aligned to the reverse strand, then the YC tag would be reverse-complemented.
 * @note YC tag that is longer than the sequence may be trimmed based on left and right hard clips extracted from the
 * CIGAR string in the BAM record.
 * @post An empty DecodedAlignments will be returned if BAM record has no YC tag or the YC tag string is empty.
 * @post Decoded alignments will be in the same orientation as the consensus sequence from the BAM record.
 * @post Left and right overhangs, if present, are assigned to R1 or R2 according to the YC tag delimiters (+/-).
 * @param record BAM record pointer
 * @return Decoded alignments
 * @throws error::Error if BAM record has an empty sequence.
 * @throws error::Error if consensus sequence and the underlying sequence represented by the YC tag still have
 * inconsistent length after trimming YC tag based on hard-clip length(s).
 * @see DecodeToBamRecords(const bam1_t* record) for decoding directly to BAM records.
 */
DecodedAlignments Decode(const bam1_t* record);

/**
 * @brief Decode a BAM record using its YC tag.
 * @note If the sequence is aligned to the reverse strand, then the YC tag would be reverse-complemented.
 * @note YC tag that is longer than the sequence may be trimmed based on left and right hard clips extracted from the
 * CIGAR string in the BAM record.
 * @post Decoded alignments will be in the same orientation as the consensus sequence from the BAM record.
 * @post Left and right overhangs, if present, are assigned to R1 or R2 according to the YC tag delimiters (+/-).
 * @post Output BAM records are generated for R1 and R2 based on the input BAM record with the following adjusted:
 * 1. read name
 * 2. reference start position
 * 3. CIGAR string
 * 4. read sequence
 * 5. base quality
 * @param record BAM record pointer
 * @return one of the following:
 * 1. empty, if no reads can be decoded
 * 2. one BAM record, if there is an overhang but no duplex segment
 * 3. two BAM records, for R1 and R2, if there is a duplex segment
 * @warning Auxiliary fields of the source BAM record are naively copied to the output BAM record(s). If the consensus
 * sequence contains discordant bases (i.e. character codes in YC tag), then auxiliary fields, such as `NM` and `MD`,
 * are very likely to be incorrect for the decoded reads. For example, a discordant code may convert a reference
 * mismatch (e.g. `NM:i:1`) from the consensus sequence into a reference match (`NM:i:0`) in a decoded read.
 * @throws error::Error if the YC tag string is malformed or cannot be parsed (e.g. invalid YC tag format or codes).
 * @throws error::Error if consensus sequence and the underlying sequence represented by the YC tag still have
 * inconsistent length after trimming YC tag based on hard-clip length(s).
 */
std::variant<std::monostate, io::Bam1Ptr, std::pair<io::Bam1Ptr, io::Bam1Ptr>> DecodeToBamRecords(const bam1_t* record);

/**
 * Ensures that the length of the read sequence matches the size of the provided base types vector.
 *
 * @param read Pointer to a BAM record (bam1_t) representing the read.
 * @param base_types A vector of yc_decode::BaseType representing the base types of the read sequence.
 *
 * @throws error::Error If the size of the base_types vector does not match the length of the read sequence.
 */
void RequireConsistentReadLength(const bam1_t* read, const vec<yc_decode::BaseType>& base_types);

/**
 * @brief Constructs an std::runtime_error indicating an invalid base type in the read with the provided read name.
 * @param read_name The name of the read where the invalid base type was encountered.
 * @param invalid_base_type The string of the invalid base type that was encountered.
 * @return std::runtime_error The constructed runtime error with the formatted message.
 */
std::runtime_error BaseTypeError(const std::string& read_name, const BaseType& invalid_base_type);

}  // namespace xoos::yc_decode
