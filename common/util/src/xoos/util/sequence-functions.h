#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>

#include <xoos/types/int.h>

namespace xoos::sequence {

/// @brief The four canonical DNA bases, used for index-to-base conversion and sequence generation.
constexpr std::string_view kDnaAlphabet = "ACGT";

/// @brief Size of the byte-indexed base lookup tables (one entry per possible unsigned char value).
constexpr std::size_t kByteTableSize = 256;

/**
 * @brief Lookup table to convert DNA/RNA base characters to a 2-bit index (A/a→0, C/c→1, G/g→2, T/t/U/u→3).
 * Accepts both uppercase and lowercase, plus uracil. All other characters map to 4 (invalid).
 */
constexpr std::array<u8, kByteTableSize> kBaseToBin = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

/**
 * @brief Lookup table to convert a DNA base to its complement, preserving case (A↔T, C↔G, a↔t, c↔g).
 * Non-ACGT characters pass through unchanged.
 */
constexpr std::array<u8, kByteTableSize> kBaseToComplement = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,
    22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,
    44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  'T',
    66,  'G', 68,  69,  70,  'C', 72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  'A', 85,  86,  87,
    88,  89,  90,  91,  92,  93,  94,  95,  96,  't', 98,  'g', 100, 101, 102, 'c', 104, 105, 106, 107, 108, 109,
    110, 111, 112, 113, 114, 115, 'a', 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131,
    132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
    154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
    198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219,
    220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241,
    242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
};

// simple reverse-complement used for processing UMIs, which may be
// set to just "*" for partial reads
std::string ReverseComplement(const std::string& str);

std::string ReverseComplementCaseSensitive(const std::string_view& str);

/// @brief Check whether a character is one of the four uppercase DNA bases (A, C, G, T).
bool IsACGT(char c);

/// @brief Check whether a character is one of the IUPAC ambiguous DNA bases (R, Y, S, W, K, M, B, D, H, V),
/// case-insensitive. Returns false for A, C, G, T, N and non-alphabetic characters.
bool IsAmbiguousBase(char c);

/**
 * Check whether an ACGT query base is consistent with an IUPAC ambiguous reference base.
 * Both arguments must be uppercase. Returns false for unknown ambiguity codes.
 * IUPAC codes: R=A|G, Y=C|T, S=G|C, W=A|T, K=G|T, M=A|C, B=C|G|T, D=A|G|T, H=A|C|T, V=A|C|G.
 */
bool AmbiguousBaseMatchesQuery(char ref_upper, char query_upper);

/**
 * Helper function to convert a base character (A, C, G, T) to its corresponding index (0, 1, 2, 3).
 * Accepts lowercase (a, c, g, t) and RNA uracil (U, u) — lowercase maps to the same index as uppercase.
 * @param base Base character (A/a, C/c, G/g, T/t/U/u)
 * @return Index (0 for A, 1 for C, 2 for G, 3 for T/U)
 * @throws std::invalid_argument if base is not a valid ACGT/U character
 */
constexpr u8 BaseToIndex(const char base) {
  const u8 bin_base = kBaseToBin.at(static_cast<u8>(base));
  if (bin_base > 3) [[unlikely]] {
    throw std::invalid_argument("Non-ACGT character in BaseToIndex");
  }
  return bin_base;
}

/**
 * Helper function to convert an index (0, 1, 2, 3) back to its corresponding base character (A, C, G, T).
 * Always returns uppercase. Round-trip with BaseToIndex is lossy for lowercase and uracil input.
 * @param index Base index (0 for A, 1 for C, 2 for G, 3 for T)
 * @return Base character corresponding to `index`
 * @throws std::invalid_argument if index is not in [0, 3]
 */
constexpr char IndexToBase(const u8 index) {
  if (index > 3) [[unlikely]] {
    throw std::invalid_argument("Index out of bounds in IndexToBase");
  }
  return kDnaAlphabet[index];
}

}  // namespace xoos::sequence
