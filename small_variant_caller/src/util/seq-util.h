#pragma once

#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include <xoos/types/int.h>

#include "xoos/types/float.h"

namespace xoos::svc {

/**
 * Synopsis:
 * This file contains utility functions for handling nucleotide sequences and variants.
 * It includes functions for converting sequences to a 2-bit representation, checking variant lengths,
 * and formatting variants.
 */

// 2-bit representation of up to 32-bp nucleotide sequence.
using TwoBit32mer = u64;
constexpr TwoBit32mer kTwoBitA{0b00};
constexpr TwoBit32mer kTwoBitC{0b01};
constexpr TwoBit32mer kTwoBitG{0b10};
constexpr TwoBit32mer kTwoBitT{0b11};

/**
 * @brief Check if a nucleotide is A, C, T, G or their lower case.
 * @param base nucleotide character
 * @return true if the check passes; false otherwise
 */
bool IsACTG(char base);

/**
 * @brief Check if a nucleotide is not A, C, T, G nor their lower case.
 * @param base nucleotide character
 * @return true if the check passes; false otherwise
 */
bool IsNotACTG(char base);

/**
 * @brief Check if a nucleotide sequence contains any character that is not A, C, T, G, nor their lower case.
 * @param seq nucleotide sequence
 * @return true if the check passes; false otherwise
 */
bool IsAnyNotACTG(std::string_view seq);

/**
 * @brief Check if a nucleotide sequence contains only A, C, T, G, or their lower case.
 * @param seq nucleotide sequence
 * @return true if the check passes; false otherwise
 */
bool ContainsOnlyACTG(std::string_view seq);

/**
 * @brief Convert a nucleotide sequence to its 2-bit representation. If the sequence is longer than 32 nt, only the
 * first 32 nt are considered.
 * @param seq nucleotide sequence
 * @return 2-bit representation of the sequence
 */
TwoBit32mer SeqToBits(const std::string& seq);

/**
 * @brief Convert a nucleotide sequence to its 2-bit representation and return it as a double.
 * @param seq nucleotide sequence
 * @return 2-bit representation of the sequence as a double
 */
f64 SeqToDouble(const std::string& seq);

/**
 * @brief Check whether a variant has a longer ALT representation than the other. Assumes both variants are at the same
 * reference position.
 * @param ref1 REF allele for variant 1
 * @param alt1 ALT allele for variant 1
 * @param ref2 REF allele for variant 2
 * @param alt2 ALT allele for variant 2
 * @return true if variant 1 has a longer ALT representation than variant 2
 */
bool HasLongerAlt(std::string_view ref1, std::string_view alt1, std::string_view ref2, std::string_view alt2);

/**
 * @brief Trim variant in a VCF record (e.g. from ATTCG->ACG to ATT->A)
 * @param ref Reference allele
 * @param alt Alternate allele
 * @return Trimmed representations of REF, ALT
 */
std::pair<std::string, std::string> TrimVariant(const std::string& ref, const std::string& alt);

/**
 * @brief Format two variants at the same position for a single VCF record
 * @param ref1 Reference allele for variant 1
 * @param alt1 Alt allele for variant 1
 * @param ref2 Reference allele for variant 2
 * @param alt2 Alt allele for variant 2
 * @return Tuple of formatted REF, ALT1, ALT2
 */
std::tuple<std::string, std::string, std::string> FormatVariants(const std::string& ref1,
                                                                 const std::string& alt1,
                                                                 const std::string& ref2,
                                                                 const std::string& alt2);

/**
 * @brief Count the number of unique k-mers in a nucleotide sequence.
 * @param seq Input nucleotide sequence
 * @param k K-mer size
 * @return Number of unique k-mers
 */
u32 CountUniqueKmers(std::string_view seq, u32 k);

/**
 * @brief Count the homopolymer length from the first base in a given nucleotide sequence.
 * A homopolymer is defined as the first base repeated at least once. Therefore, the output is either 0 or >= 2.
 * @param seq Input nucleotide sequence
 * @return Homopolymer length
 */
u32 GetHomopolymerLength(std::string_view seq);

/**
 * @brief Count consecutive repeats of the first k-mer in a given nucleotide sequence.
 * A repeat is defined as the first k-mer repeated at least once. Therefore, the output is either 0 or > 1.
 * @param seq Input nucleotide sequence
 * @param k K-mer size
 * @return Number of occurrences of the repeating k-mer
 */
u32 CountRepeats(std::string_view seq, u32 k);

/**
 * @brief Find the start position of the homopolymer overlapping the current position.
 * @param seq nucleotide sequence
 * @param pos current position
 * @param min_hp_length minimum length of homopolymer to be considered
 * @param min_start minimum start position of homopolymer
 * @return start position of homopolymer if present, std::nullopt otherwise
 * @note min_hp_length >= 2 because a homopolymer requires at least two identical bases
 */
std::optional<u64> FindOverlappingHomopolymer(const std::string& seq, u64 pos, u64 min_hp_length, u64 min_start);

/**
 * @brief Calculate the GC content (fraction of G and C bases or their lower case) in a nucleotide sequence.
 * @param seq nucleotide sequence
 * @return GC content as a fraction in [0, 1], or 0 if the sequence is empty or contains no ACTG characters
 * @note Non-ACTG characters are excluded from both the numerator and denominator.
 */
f32 CalculateGcContent(std::string_view seq);

}  // namespace xoos::svc
