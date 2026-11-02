#include "rdb-2bit-utils.h"

#include <algorithm>
#include <string>
#include <vector>

#include <xoos/util/sequence-functions.h>

#include "sequence/encoding-constants.h"

namespace xoos::demux {

constexpr std::vector<u8> GenerateReverseLookupTable() {
  std::vector<u8> lookup_table(static_cast<size_t>(kByteValueCount));

  // Reverse the order of the four 2-bit bases packed into each byte, i.e. base positions 0<->3 and 1<->2.
  // The masks select each base position and the shifts move it to its mirrored position.
  constexpr u32 kBase0Mask = 0b0000'0011U;
  constexpr u32 kBase1Mask = 0b0000'1100U;
  constexpr u32 kBase2Mask = 0b0011'0000U;
  constexpr u32 kBase3Mask = 0b1100'0000U;
  constexpr u32 kBase0ToBase3Shift = 6U;
  constexpr u32 kBase1ToBase2Shift = 2U;
  constexpr u32 kBase2ToBase1Shift = 2U;
  constexpr u32 kBase3ToBase0Shift = 6U;

  for (u32 i = 0; i < kByteValueCount; ++i) {
    lookup_table[i] = static_cast<u8>((i & kBase0Mask) << kBase0ToBase3Shift | (i & kBase1Mask) << kBase1ToBase2Shift |
                                      (i & kBase2Mask) >> kBase2ToBase1Shift | (i & kBase3Mask) >> kBase3ToBase0Shift);
  }

  return lookup_table;
}

static const std::vector<u8> kReverseLookupTable = GenerateReverseLookupTable();

void Reverse2BitOrder(u8* sequence, std::size_t sequence_length) {
  for (std::size_t i = 0; i < sequence_length; ++i) {
    sequence[i] = kReverseLookupTable[sequence[i]];
  }
}

constexpr std::vector<std::string> GenerateDnaDecodeTable() {
  std::vector<std::string> lookup_table(static_cast<size_t>(kByteValueCount));
  for (std::size_t i = 0; i < lookup_table.size(); ++i) {
    std::string decoded_bases;
    for (u32 j = 0; j < kBasesPerByte; ++j) {
      decoded_bases += sequence::kDnaAlphabet[(i >> (kBitsPerEncodedBase * j)) & kBaseMask];
    }
    lookup_table[i] = decoded_bases;
  }
  return lookup_table;
}

static const std::vector<std::string> kDnaLookupTable = GenerateDnaDecodeTable();

void DecodeDnaBases(const u8* const encoded, char* const buffer, const std::size_t length) {
  char* out = buffer;
  for (std::size_t i = 0; i < length >> kBitsPerEncodedBase; ++i) {
    const auto byte = encoded[i];
    const std::string& decoded_bases = kDnaLookupTable[byte];
    out = std::copy_n(decoded_bases.c_str(), decoded_bases.size(), out);
  }
  const std::size_t remaining = length & kBaseMask;
  if (remaining != 0) {
    const auto byte = encoded[length >> kBitsPerEncodedBase];
    const std::string& decoded_bases = kDnaLookupTable[byte];
    std::copy_n(decoded_bases.c_str(), remaining, out);
  }
}

constexpr std::vector<std::string> GenerateQualDecodeTable() {
  constexpr char kQual[] = {'+', '5', '?', 'I'};

  std::vector<std::string> lookup_table(static_cast<size_t>(kByteValueCount));
  for (u32 i = 0; i < kByteValueCount; ++i) {
    std::string decoded_quals;
    for (s32 j = kBasesPerByte - 1; j >= 0; --j) {
      decoded_quals += kQual[(i >> (kBitsPerEncodedBase * static_cast<u32>(j))) & kBaseMask];
    }
    lookup_table[i] = decoded_quals;
  }
  return lookup_table;
}

static const std::vector<std::string> kQualLookupTable = GenerateQualDecodeTable();

void DecodeQualScores(const u8* const encoded, char* const buffer, const std::size_t length) {
  char* out = buffer;
  for (std::size_t i = 0; i < length >> kBitsPerEncodedBase; ++i) {
    uint8_t byte = encoded[i];
    const std::string& decoded_qual = kQualLookupTable[byte];
    out = std::copy_n(decoded_qual.c_str(), decoded_qual.size(), out);
  }
  const auto remaining = length & kBaseMask;
  if (remaining != 0) {
    const auto byte = encoded[length >> kBitsPerEncodedBase];
    const std::string& decoded_qual = kQualLookupTable[byte];
    std::copy_n(decoded_qual.c_str(), remaining, out);
  }
}

}  // namespace xoos::demux
