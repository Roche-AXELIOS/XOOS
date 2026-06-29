/**
 * This file (rdb-2bit-utils.h) contains utility functions to decode 2-bit
 * encoded DNA bases and quality scores.
 */

#pragma once

#include <xoos/types/int.h>

namespace xoos::demux {
void Reverse2BitOrder(u8* sequence, std::size_t sequence_length);
void DecodeDnaBases(const u8* encoded, char* buffer, std::size_t length);
void DecodeQualScores(const u8* encoded, char* buffer, std::size_t length);
}  // namespace xoos::demux
