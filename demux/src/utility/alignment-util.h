#pragma once

#include <xoos/types/int.h>

#include <string_view>

namespace xoos::demux {

// quality base types
constexpr char kDiscordBaseQual{'&'}, kSimplexBaseQual{'7'}, kConcordBaseQual{'H'};

// Edit operations.
enum { kEdopMatch = 0, kEdopInsert = 1, kEdopDelete = 2, kEdopMismatch = 3, kEdopMatchMethyl = 4 };

// alphabet for DNA bases conversion
constexpr char kAlphabet[4] = {'A', 'C', 'G', 'T'};

// YC tag components
constexpr char kDefaultJunction{'+'};
constexpr char kSwappedJunction{'-'};

// IUPAC-like base encoding table for duplex consensus. This table is used to encode two bases into one IUPAC-like base
// and is used to create a compact, lossless compression of regions R1 and R2 used in the alignment process. The
// IUPAC-like encoding describes the differences between each of the two reads and the consensus sequence and is
// therefore a 5x5 table because we also need to encode indels.
//                                    A    C    G    T    -
constexpr char kIUPACTable[5][5] = {{'*', 'M', 'R', 'W', 'I'},   // A
                                    {'B', '*', 'S', 'Y', 'P'},   // C
                                    {'D', 'V', '*', 'K', 'J'},   // G
                                    {'H', 'E', 'F', '*', 'X'},   // T
                                    {'L', 'Q', 'O', 'Z', '*'}};  // -

// Index in kIUPACTable for indels, which is the last row/column
constexpr u8 kIndelIndex = 4;

void LeftAlignConvert(std::string_view consensus, u8* alignment);

void GroupRightAlignConvert(std::string_view consensus, u8* alignment);

void MethylConvert(char* consensus, u8* alignment, s32 alignment_length, char* xm_tag_start, const u8* r2,
                   bool swapped);

};  // namespace xoos::demux
