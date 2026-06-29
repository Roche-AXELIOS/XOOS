#include "rdb-2bit-utils.h"

#include <cstring>
#include <string>
#include <vector>

namespace xoos::demux {

constexpr std::vector<u8> GenerateReverseLookupTable() {
  std::vector<u8> lookup_table(256);

  for (int i = 0; i < 256; ++i) {
    lookup_table[i] = (i & 0b11) << 6 | (i & 0b1100) << 2 | (i & 0b110000) >> 2 | (i & 0b11000000) >> 6;
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
  constexpr char kDNA[] = {'A', 'C', 'G', 'T'};

  std::vector<std::string> lookup_table(256);
  for (std::size_t i = 0; i < lookup_table.size(); ++i) {
    std::string decoded_bases;
    for (int j = 0; j < 4; ++j) {
      decoded_bases += kDNA[(i >> (2 * j)) & 0b11];
    }
    lookup_table[i] = decoded_bases;
  }
  return lookup_table;
}

static const std::vector<std::string> kDnaLookupTable = GenerateDnaDecodeTable();

void DecodeDnaBases(const u8* encoded, char* buffer, std::size_t length) {
  char* out = buffer;
  for (std::size_t i = 0; i < length >> 2; ++i) {
    uint8_t byte = encoded[i];
    const std::string& decoded_bases = kDnaLookupTable[byte];
    std::memcpy(out, decoded_bases.c_str(), decoded_bases.size());
    out += decoded_bases.size();
  }
  std::size_t remaining = length & 0b11;
  if (remaining != 0) {
    uint8_t byte = encoded[length >> 2];
    const std::string& decoded_bases = kDnaLookupTable[byte];
    std::memcpy(out, decoded_bases.c_str(), remaining);
    out += remaining;
  }
}

constexpr std::vector<std::string> GenerateQualDecodeTable() {
  constexpr char kQual[] = {'+', '5', '?', 'I'};

  std::vector<std::string> lookup_table(256);
  for (int i = 0; i < 256; ++i) {
    std::string decoded_quals;
    for (int j = 3; j >= 0; --j) {
      decoded_quals += kQual[(i >> (2 * j)) & 0b11];
    }
    lookup_table[i] = decoded_quals;
  }
  return lookup_table;
}

static const std::vector<std::string> kQualLookupTable = GenerateQualDecodeTable();

void DecodeQualScores(const u8* encoded, char* buffer, std::size_t length) {
  char* out = buffer;
  for (std::size_t i = 0; i < length >> 2; ++i) {
    uint8_t byte = encoded[i];
    const std::string& decoded_qual = kQualLookupTable[byte];
    std::memcpy(out, decoded_qual.c_str(), decoded_qual.size());
    out += decoded_qual.size();
  }
  std::size_t remaining = length & 0b11;
  if (remaining != 0) {
    uint8_t byte = encoded[length >> 2];
    const std::string& decoded_qual = kQualLookupTable[byte];
    std::memcpy(out, decoded_qual.c_str(), remaining);
    out += remaining;
  }
}

}  // namespace xoos::demux
