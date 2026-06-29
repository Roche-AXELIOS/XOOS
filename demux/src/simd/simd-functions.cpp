#include "simd/simd-functions.h"

#include <xoos/util/char-functions.h>

#include <algorithm>
#include <utility>

// Generates code for every target that this compiler can support. The foreach_target.h
// file will cause this source to be compiled multiple times for each target; we need to
// pass the name of this file to the compiler so that it can be included in the preprocessor
// - hence the inclusion of the full path of this file (as seen by the Highway library).
#undef HWY_TARGET_INCLUDE
#define HWY_DISABLED_TARGETS (HWY_SVE | HWY_SVE2 | HWY_SVE_256 | HWY_SVE2_128)
// Defined in CMakeLists.txt file so build will be consistent from within Docker / CI or natively.
#define HWY_TARGET_INCLUDE SIMD_FUNCTIONS_FILE

#include <hwy/foreach_target.h>  // must come before highway.h
#include <hwy/highway.h>
#include <hwy/print-inl.h>

#include "adapters/duplex/duplex-match.h"
#include "io/read-record.h"

// Required for ARM64
using hwy::HWY_NAMESPACE::And;
using hwy::HWY_NAMESPACE::TableLookupLanes;

namespace HWY_NAMESPACE {  // required: unique per target

constexpr HWY_FULL(uint64_t) kFullWidth64;
constexpr HWY_FULL(uint32_t) kFullWidth32;
constexpr HWY_FULL(uint16_t) kFullWidth16;
constexpr HWY_FULL(uint8_t) kFullWidth8;
constexpr HWY_FULL(int64_t) kFullWidth64S;
constexpr HWY_FULL(int32_t) kFullWidth32S;

// The methodology used for 2-bit encoding is slightly different if we want to generate complementary bits; the
// target encoding is A=00, C=01, G=10, T=11. We would like to determine the encoding using the ascii values of the
// input, which are: 'A': 65, 'C': 67, 'G': 71, 'T': 84, which is
// A: 01000001
// C: 01000011
// G: 01000111
// T: 01010100
// We can see that bit 2 provides us one of the output bits (the most significant bit). The least significant bit is
// trickier, but staring at the values, we can see that the least significant bit can be obtained by XOR'ing bit 1 and
// bit 2 of the input.

HWY_ATTR void To2BitSIMD(const uint8_t* p_seq, uint32_t length, uint8_t* p_out) {
  constexpr HWY_FULL(uint8_t) kFullWidth;  // Specify the maximum width to process; this will be 64 bytes in AVX3, but
                                           // we only will load 32 bytes (2 bytes for every base during processing).

  // After reading the input sequence, we need to replicate each base two times; we'll use an index into a lookup
  // table for that. This maps very well on AVX3 hardware.
  alignas(64) constexpr uint8_t kIndices[] = {0,  0,  1,  1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,
                                              8,  8,  9,  9,  10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15,
                                              16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23,
                                              24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31};

  // Calculate the number of bases that can be processed in one go.
  const auto incr{static_cast<int>(Lanes(kFullWidth) >> 1)};  // this will be 32 in AVX3
  const auto incr2{incr >> 2};  // increment in output (2 bits per base, 4 bases per byte, is shift by 2)
  const auto indices = SetTableIndices(kFullWidth, kIndices);  // Load the indices to be loaded into a register

  constexpr HWY_CAPPED(uint16_t, 32)
      kCappedInt;  // specify a register for 16-bit integers; cap max lanes at 32 (half of max AVX3)
  const auto mask1{Set(kFullWidth, 0x04)};  // Extract bit #2 from the input
  const auto mask2{
      ResizeBitCast(kFullWidth, Set(kCappedInt, 0x0004))};  // Extract bit #2 from the input (every other input)

  // Split the processing into a loop that processes as many bases as possible (up to 32 bases at a time),
  // followed by a wrap-up part that processes the remaining bases.
  int offset{0}, out_offset{0};
  const auto stop_offset{static_cast<int>(length) - 2 * incr};
  for (; offset <= stop_offset; offset += incr, out_offset += incr2) {  // process the maximum possible number of bases
    const auto in = LoadU(kFullWidth, p_seq + offset);                  // load up to 32 bytes
    const auto in2 = TableLookupLanes(in, indices);                     // Duplicate every base.
    const auto m1 = And(mask1, in2);                                    // Extract bit 2
    const auto in3 = ShiftLeft<1>(in2);  // Shift data left by 1 to get bit 1 into the position of bit 2
    const auto m2 = And(mask2, in3);     // Extract bit 1 at position of bit 2
    const auto result{Xor(m1, m2)};      // XOR will give us the bit pattern we want (only bit 2 is relevant)
    // AND the sequence with the mask and compare with the mask, then write out the mask
    StoreMaskBits(kFullWidth, Eq(result, mask1), p_out + out_offset);  // Write the resulting mask
  }

  // Tail handling: last (incr + (length%incr)) bases are processed in two iterations.
  // iteration 1: incr bases
  // iteration 2: last length % incr bases
  // iteration 1 could not be included in the above loop as that would cause out of bounds access.
  while (std::cmp_less(offset, length)) {
    const auto num_to_load = std::min(static_cast<int>(length - offset), incr);
    const auto in = LoadN(kFullWidth, p_seq + offset, num_to_load);  // load remaining bases
    const auto in2 = TableLookupLanes(in, indices);                  // Duplicate every base.
    const auto m1 = And(mask1, in2);                                 // Extract bit 2
    const auto in3 = ShiftLeft<1>(in2);  // Shift left by 1 to get bit 1 into the position of bit 2
    const auto m2 = And(mask2, in3);     // Extract bit 1
    const auto result{Xor(m1, m2)};      // XOR will give us the bit pattern we want (only bit 2 is relevant)

    // Assume that we do not have to worry about access violations when writing data; the caller should have
    // added enough padding to the output buffer that we do not worry about it.
    StoreMaskBits(kFullWidth, Eq(result, mask1), p_out + out_offset);  // Write the resulting mask
    offset += num_to_load;
    out_offset += (num_to_load >> 2);
  }
}

HWY_ATTR uint64_t FindSpaceSIMD(const uint8_t* const buf, const uint64_t begin, const uint64_t end) {
  constexpr HWY_FULL(uint8_t) kFullWidth;
  const auto incr{Lanes(kFullWidth)};
  const auto mask_sp = Set(kFullWidth, ' ');    // Load the space character into a register
  const auto mask_tab = Set(kFullWidth, '\t');  // Load the tab character into a register
  const auto mask_cr = Set(kFullWidth, '\r');   // Load the carriage return character into a register

  auto i = begin;
  for (; i + incr <= end; i += incr) {
    const auto in = LoadU(kFullWidth, buf + i);  // Load up to 64 bytes
    // Check for space, tab, line feed, carriage return, vertical tab, and form feed
    const auto mask = Or(Eq(in, mask_sp), And(Ge(in, mask_tab), Le(in, mask_cr)));

    if (!AllFalse(kFullWidth, mask)) {
      return i + FindKnownFirstTrue(kFullWidth, mask);
    }
  }
  // If we arrive here, no space was found in the last bytes (usually 32 or 64). We need to check the remaining bytes.
  for (; i < end; ++i) {
    if (xoos::character::IsSpace(buf[i])) {
      return i;
    }
  }
  return end;
}

HWY_ATTR uint64_t FindNonTrivialSpaceSIMD(const uint8_t* const buf, const uint64_t begin, const uint64_t end) {
  constexpr HWY_FULL(uint8_t) kFullWidth;
  const auto incr{Lanes(kFullWidth)};
  const auto mask_tab = Set(kFullWidth, '\t');  // Load the tab character into a register
  const auto mask_cr = Set(kFullWidth, '\r');   // Load the carriage return character into a register

  auto i = begin;
  for (; i + incr <= end; i += incr) {
    const auto in = LoadU(kFullWidth, buf + i);  // Load up to 64 bytes
    // Check for tab, line feed, carriage return, vertical tab, and form feed
    const auto mask = And(Ge(in, mask_tab), Le(in, mask_cr));

    if (!AllFalse(kFullWidth, mask)) {
      return i + FindKnownFirstTrue(kFullWidth, mask);
    }
  }
  // If we arrive here, no space was found in the last bytes (usually 32 or 64). We need to check the remaining bytes.
  for (; i < end; ++i) {
    if (xoos::character::IsNonTrivialSpace(buf[i])) {
      return i;
    }
  }
  return end;
}

// The edlib library assumes a kWordSize of 64 (i.e. 64 bits per register).
constexpr int kWordSize = 64;

inline int CeilDiv64(int x) { return (x + kWordSize - 1) >> 6; }  // divide by 64 = shift by 6

// This code is used for pairwise alignment; it is an optimized version of the TransformSequences()
// function that are part of the edlib library, of which an optimized version is embedded in our Demux
// library.
HWY_ATTR void TransformSequenceSIMD(const unsigned char* buf, int begin, int end, uint8_t* digit, bool reverse,
                                    bool invert) {
  // First step: convert the input sequence into a 8-bit encoding (A=0, C=1, G=2, T=3). Instead of a LUT, we'll use
  // the bit pattern of the ASCII values to determine the encoding. This is a very efficient way to encode the data.
  // We have two code paths: one for the forward sequence and one for the reverse sequence.
  // 'A': binary 01000001     -> needs to map on 00
  // 'C': binary 01000011     -> need to map on 01
  // 'G': binary 01000111     -> needs to map on 10
  // 'T': binary 01010100     => needs to map on 11
  // We can see that bit 2 provides us one of the output bits (the most significant bit). The least significant bit is
  // trickier, but staring at the values, we can see that the least significant bit can be obtained by XOR'ing bit 1 and
  // bit 2 of the input.
  // Calculate the number of bases that can be processed in one go.
  const auto incr{static_cast<int>(Lanes(kFullWidth8))};       // this will be 32 in AVX3
  const auto mask_not_bit0{Set(kFullWidth8, 0x01)};            // Zero all bits except bit #0
  const auto mask_not_bit1{Set(kFullWidth8, 0x02)};            // Zero all bits except bit #1
  const auto invert_mask{Set(kFullWidth8, invert ? 0x3 : 0)};  // Invert mask if needed

  // We loop over the OUTPUT that we need to generate.
  int number_bases{1 + end - begin};

  // First, convert the input into a 8-bit encoding with 4 possible values (0..3).
  int offset;
  for (offset = 0; offset + incr <= number_bases; offset += incr) {
    auto in{Undefined(kFullWidth8)};  // will contain the input data.
    if (reverse) {
      in = LoadU(kFullWidth8, buf + end - incr + 1 - offset);  // last byte should become the first one
      in = Reverse(kFullWidth8, in);                           // reverse the order of input
    } else {
      in = LoadU(kFullWidth8, buf + begin + offset);  // load up to 32 bytes
    }

    const auto shift1{ShiftRight<1>(in)};         // Shift data right by 1 to get bit 1 into the position of bit 0
    const auto shift2{ShiftRight<2>(in)};         // Shift data right by 2 to get bit 2 into the position of bit 0
    const auto xor0{Xor(shift1, shift2)};         // Xor of bit1 and bit2 is least-significant bit of 2-bit encoding
    const auto bit1{And(mask_not_bit1, shift1)};  // Masking out bit 2 of input -> bit 1 of output
    const auto bit0{And(mask_not_bit0, xor0)};    // Masking out bit 1 of xor'd bits -> bit 0 of output
    const auto digit_out{Or(bit0, bit1)};         // Combine the two bits to get the 2-bit encoding
    StoreU(Xor(digit_out, invert_mask), kFullWidth8, digit + offset);  // Store the (inverted) 8-bit encoding
  }

  auto num_bases_remaining = number_bases - offset;
  if (num_bases_remaining > 0) {
    auto in{Undefined(kFullWidth8)};
    if (reverse) {
      in = LoadN(kFullWidth8, buf + begin, num_bases_remaining);
      in = SlideUpLanes(kFullWidth8, in, incr - num_bases_remaining);
      in = Reverse(kFullWidth8, in);
    } else {
      in = LoadN(kFullWidth8, buf + begin + offset, num_bases_remaining);
    }

    const auto shift1{ShiftRight<1>(in)};         // Shift data right by 1 to get bit 1 into the position of bit 0
    const auto shift2{ShiftRight<2>(in)};         // Shift data right by 2 to get bit 2 into the position of bit 0
    const auto xor0{Xor(shift1, shift2)};         // Xor of bit1 and bit2 is least-significant bit of 2-bit encoding
    const auto bit1{And(mask_not_bit1, shift1)};  // Masking out bit 2 of input -> bit 1 of output
    const auto bit0{And(mask_not_bit0, xor0)};    // Masking out bit 1 of xor'd bits -> bit 0 of output
    const auto digit_out{Or(bit0, bit1)};         // Combine the two bits to get the 2-bit encoding
    StoreN(Xor(digit_out, invert_mask), kFullWidth8, digit + offset,
           num_bases_remaining);  // Store the (inverted) 8-bit encoding
  }
}

HWY_ATTR void ReverseCopySIMD(const unsigned char* p_start, const unsigned char* p_end, unsigned char* p_out) {
  // Optimized version of std::reverse_copy() which is not fast enough for our purposes.
  const auto incr{static_cast<int>(Lanes(kFullWidth8))};  // this will be 64 in AVX3
  auto length{static_cast<int>(p_end - p_start)};

  int i = 0;
  for (; i < (length - incr); i += incr) {
    const auto in{LoadU(kFullWidth8, p_end - i - incr)};
    const auto reversed{Reverse(kFullWidth8, in)};
    StoreU(reversed, kFullWidth8, p_out + i);
  }
  // Remaining bytes need to be copied one by one.
  for (; i < length; ++i) {
    p_out[i] = p_end[-i - 1];
  }
}

// Same code as above, but now the Transform step is omitted and only BuildPeq() is being done
HWY_ATTR void BuildPeqSIMD(const uint8_t* digit, int number_bases, uint64_t* bit) {
  const auto incr{static_cast<int>(Lanes(kFullWidth8))};  // this will be 64 in AVX3

  // Now create the PEQ table. Instead of creating this table "on-the-fly" while converting the input to an 8-bit
  // representation, we will loop over that 8-bit representation multiple times for each binary output we need to
  // create. This is a speed optimization; as we need to write multiple outputs where we only calculate 32 bits at
  // a time (4 bytes), writing to multiple outputs will be slow because we can't take advantage of write-combining.
  // It's faster to write to a single output at a time, the loop overhead is relatively small.
  // The PEQ table contains bit masks for every base type (A, C, G, T). While it is passed as a 64-bit uint64_t array,
  // it is actually an uint64_t[5][dynamic] array where any unused entries are set to 1.
  int max_num_blocks = CeilDiv64(number_bases);
  const auto incr2{incr >> 3};  // 8 bases per byte
  for (int base = 0; base < 4; ++base) {
    const auto cmp_base{Set(kFullWidth8, base)};
    const auto output{reinterpret_cast<uint8_t*>(bit + base * max_num_blocks)};
    auto offset2{0};

    for (int offset = 0; offset < number_bases; offset += incr, offset2 += incr2) {
      const auto in{LoadU(kFullWidth8, digit + offset)};   // load up to 32 bytes
      const auto mask{Eq(in, cmp_base)};                   // Compare the input with the base
      StoreMaskBits(kFullWidth8, mask, output + offset2);  // Store the PEQ table for the base
    }

    // Remaining bases need all to be set to ones.
    auto remaining{number_bases & (kWordSize - 1)};  // number of bases that are not a multiple of 64
    if (remaining)                                   // if number_bases is not a multiple of 64
    {
      // If remaining is 1, we need to set the last 63 bits to 1.
      // If remaining is 2, we need to set the last 62 bits to 1, etc.
      uint64_t mask{(1ULL << remaining) - 1};
      bit[base * max_num_blocks + max_num_blocks - 1] |= ~mask;
    }
  }

  // We're almost done. The last step is setting the PEQ table for the "N" base - which is actually not used, I
  // believe (might be omitted soon).
  for (int i = 0; i < max_num_blocks; ++i) {
    bit[4 * max_num_blocks + i] = ~0ULL;  // Set all bits to 0
  }
}

constexpr uint64_t kReverseBases[]{
    0x0,  0x40, 0x80, 0xc0, 0x10, 0x50, 0x90, 0xd0, 0x20, 0x60, 0xa0, 0xe0, 0x30, 0x70, 0xb0, 0xf0, 0x4,  0x44, 0x84,
    0xc4, 0x14, 0x54, 0x94, 0xd4, 0x24, 0x64, 0xa4, 0xe4, 0x34, 0x74, 0xb4, 0xf4, 0x8,  0x48, 0x88, 0xc8, 0x18, 0x58,
    0x98, 0xd8, 0x28, 0x68, 0xa8, 0xe8, 0x38, 0x78, 0xb8, 0xf8, 0xc,  0x4c, 0x8c, 0xcc, 0x1c, 0x5c, 0x9c, 0xdc, 0x2c,
    0x6c, 0xac, 0xec, 0x3c, 0x7c, 0xbc, 0xfc, 0x1,  0x41, 0x81, 0xc1, 0x11, 0x51, 0x91, 0xd1, 0x21, 0x61, 0xa1, 0xe1,
    0x31, 0x71, 0xb1, 0xf1, 0x5,  0x45, 0x85, 0xc5, 0x15, 0x55, 0x95, 0xd5, 0x25, 0x65, 0xa5, 0xe5, 0x35, 0x75, 0xb5,
    0xf5, 0x9,  0x49, 0x89, 0xc9, 0x19, 0x59, 0x99, 0xd9, 0x29, 0x69, 0xa9, 0xe9, 0x39, 0x79, 0xb9, 0xf9, 0xd,  0x4d,
    0x8d, 0xcd, 0x1d, 0x5d, 0x9d, 0xdd, 0x2d, 0x6d, 0xad, 0xed, 0x3d, 0x7d, 0xbd, 0xfd, 0x2,  0x42, 0x82, 0xc2, 0x12,
    0x52, 0x92, 0xd2, 0x22, 0x62, 0xa2, 0xe2, 0x32, 0x72, 0xb2, 0xf2, 0x6,  0x46, 0x86, 0xc6, 0x16, 0x56, 0x96, 0xd6,
    0x26, 0x66, 0xa6, 0xe6, 0x36, 0x76, 0xb6, 0xf6, 0xa,  0x4a, 0x8a, 0xca, 0x1a, 0x5a, 0x9a, 0xda, 0x2a, 0x6a, 0xaa,
    0xea, 0x3a, 0x7a, 0xba, 0xfa, 0xe,  0x4e, 0x8e, 0xce, 0x1e, 0x5e, 0x9e, 0xde, 0x2e, 0x6e, 0xae, 0xee, 0x3e, 0x7e,
    0xbe, 0xfe, 0x3,  0x43, 0x83, 0xc3, 0x13, 0x53, 0x93, 0xd3, 0x23, 0x63, 0xa3, 0xe3, 0x33, 0x73, 0xb3, 0xf3, 0x7,
    0x47, 0x87, 0xc7, 0x17, 0x57, 0x97, 0xd7, 0x27, 0x67, 0xa7, 0xe7, 0x37, 0x77, 0xb7, 0xf7, 0xb,  0x4b, 0x8b, 0xcb,
    0x1b, 0x5b, 0x9b, 0xdb, 0x2b, 0x6b, 0xab, 0xeb, 0x3b, 0x7b, 0xbb, 0xfb, 0xf,  0x4f, 0x8f, 0xcf, 0x1f, 0x5f, 0x9f,
    0xdf, 0x2f, 0x6f, 0xaf, 0xef, 0x3f, 0x7f, 0xbf, 0xff};

inline uint64_t LoadReverse64(uint8_t const* b) {
  union {
    uint64_t u64;
    uint8_t u8[8];
  } u = {.u8 = {b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]}};

  uint64_t out = kReverseBases[u.u8[0]];
  out = (out << 8) | kReverseBases[u.u8[1]];
  out = (out << 8) | kReverseBases[u.u8[2]];
  out = (out << 8) | kReverseBases[u.u8[3]];
  out = (out << 8) | kReverseBases[u.u8[4]];
  out = (out << 8) | kReverseBases[u.u8[5]];
  out = (out << 8) | kReverseBases[u.u8[6]];
  out = (out << 8) | kReverseBases[u.u8[7]];

  return out;
}

HWY_ATTR int FindSymmetryPositionSIMD(xoos::demux::FixedReadRecord& record) {
  // Step 1: Find a point of interest in the 3p end of the read. We're using 2-bit encoding, so let's start 32 bases
  // from the end of the read. This is a good distance to avoid any junk that might be present at the end of the read.
  // So we'll first evaluate bases 32-47 from the end, then 40-55, then 48-63, etc.

  // Recent Xeon processors introduced support for affine Galois field transformations, which can be used to reverse
  // bits efficiently. To trigger that code path in Highway, we need to load the data into a vectorized type, it does
  // not really matter which type we use as we only will use 32 bits of it. Note: position needs to be a multiple
  // of 4.

  // Reverse_in_3p now has 64 bits reversed, but we only need the 32 highest bits of it. We need to extract and
  // duplicate the "even" 32-bits of the register.
  // Do one iteration to start with
  size_t pos_5p = 0;  // 8 bases from the 5p start, don't need the preamble

  auto length{record.SeqLen()};
  uint32_t min_pos{(length * 7) >> 4};  // 5p needs to be larger than 7/16 of the length of the read
  const uint8_t* start_ptr{record.TwoBitsSeq()};

  // Create a threshold mask
  const auto threshold_mask{Set(kFullWidth32, ~0)};
  const size_t max_search{length >> 2};  // don't forget to divide by 4

  // This is not an exhaustive search, but a "most bang for the buck" search. If we don't find symmetry,
  // we'll execute an alternate strategy so we can afford to miss some matches.
  for (; pos_5p < max_search; pos_5p += Lanes(kFullWidth8)) {
    auto in5p = LoadU(kFullWidth8, start_ptr + pos_5p);  // load up to 32 bytes (128 bases)

    // This is a sliding window search. We can affect the results by changing the search range max_search and the
    // increment in the offset. No improvement by using an increment smaller than 4.
    constexpr auto kIncrement{4};  // sweet spot

    size_t start_search{(length - 32) >> 2};
    for (size_t o = 0; o < max_search && o < start_search; o += kIncrement) {
      // Load a 64-bit value (32 bases) from position pos_3p + o, then reverse the bases. This causes the first base
      // to be located at the highest 2 bits of the 64-bit value.
#if HWY_TARGET <= HWY_AVX3_DL
      auto in128{LoadU(kFullWidth8, start_ptr + start_search - o)};               // SIMD load 128-bit value
      auto in64{Broadcast<0>(hwy::HWY_NAMESPACE::BitCast(kFullWidth64, in128))};  // SIMD load 64-bit value
      //  auto in64 { Set(kFullWidth64,simd::Load64(start_ptr + start_search - o)) };   // SIMD load 64-bit value
      auto in64reverse = hwy::HWY_NAMESPACE::ReverseLaneBytes(hwy::HWY_NAMESPACE::BitCast(
          kFullWidth64,
          hwy::HWY_NAMESPACE::detail::GaloisAffine(hwy::HWY_NAMESPACE::BitCast(kFullWidth8, in64),
                                                   hwy::HWY_NAMESPACE::Set(kFullWidth64, 0x4080102004080102u))));
      // After reversing the 64-bit value, we need to extract the 32 highest bits of it.

      for (int offset = 0; offset < 16; ++offset) {
        const auto in32 = BroadcastLane<1>(BitCast(kFullWidth32, in64reverse));
        in64reverse = ShiftLeft<2>(in64reverse);  // shift the data to the left by 1 base (2 bits)
#else
      auto reverse_64{LoadReverse64(start_ptr + start_search - o)};

      for (int offset = 0; offset < 16; ++offset) {
        union {
          uint64_t u64;
          uint32_t u32[2];
        } u = {.u64 = reverse_64};

        auto reverse32 =
            u.u32[1];  // get the highest 32 bits, so these 16 bases are associated with the first position.
        const auto in32{Set(kFullWidth32, reverse32)};
        reverse_64 <<= 2;  // shift the data to the left by 1 base (2 bits)
#endif

        // Here's the magic: perform an XOR of the 32-bit data with the 32-bit data from the 5p end of the read. For
        // a complementary match, the XOR will give us a 32-bit value of 0xFFFFFFFF.
        const auto in1{BitCast(kFullWidth32, in32)};
        const auto in2{BitCast(kFullWidth32, in5p)};
        auto xor_result{Xor(in1, in2)};

        auto hits{Eq(xor_result, threshold_mask)};

        while (!AllFalse(kFullWidth32, hits)) {
          // Find the first uint32 register that has a complete match, expect value is [0..15]. Each register contains
          // 16 bases, so we need to multiply the value by 16 (shift left by 4) to get the position of the match.

          auto index = FindKnownFirstTrue(kFullWidth32, hits);
          auto pf = index << 4;  // first int register with a complete match

          // We have a match, now we need to find the position of the match.
          // On the 5p side, we read data at position pos_5p, so as every byte contains 4 bases, we need to
          // multiply the position by 4. The variable p gives us the first int32 register that contained a complete
          // match, every register contains 16 bases. This allows us to calculate the 5p position of the match:
          record.trim_info_duplex.symmetry_5p_pos = (pos_5p << 2) + pf;  // position converted back to bases

          // At the 3p end, we read 64-bit (32 bases) of data at position start_search - o. I determined the formula
          // for the position of the match by trial and error. The formula is: ((start_search - o) << 2) + 15 + offset.
          record.trim_info_duplex.symmetry_3p_pos = ((start_search - o) << 2) + 15 + offset;

          record.trim_info_duplex.symmetry_middle =
              (record.trim_info_duplex.symmetry_5p_pos + record.trim_info_duplex.symmetry_3p_pos) >> 1;

          if (record.trim_info_duplex.symmetry_3p_pos < length && record.trim_info_duplex.symmetry_5p_pos < length &&
              record.trim_info_duplex.symmetry_middle > min_pos) {
            return static_cast<int>(record.trim_info_duplex.symmetry_middle);
          } else {  // this position was not usable (symmetry falls outside the read, probably caused by spurious
                    // remnants
            // of previous reads). We'll try again with a different position.
            // pf is the index position of the "hit". We need to set that value to zero now. By doing this, we'll
            // catch a potential match at the next position (the number of reads affected is very low and one could
            // just return '-1' instead without noticeable impact - but throughput is basically identical so let's do
            // it 'right').
            xor_result = InsertLane(xor_result, index, 0u);
            hits = Eq(xor_result, threshold_mask);
          }
        }
      }
    }
  }

  return -1;
}

// Finds the leading run of zero-valued bytes (edlib match encoding) using SIMD comparison,
// then copies the corresponding bases from the input sequence to the consensus buffer.
// Uses trailing zero count on the inverted match mask to determine the run length efficiently.
HWY_ATTR int CopyMatchingSeqSIMD(const uint8_t* alignment, uint8_t* bases_in, uint8_t* consensus_out) {
  size_t nr_bases{0};
  const auto incr{Lanes(kFullWidth8)};  // this will be 64 in AVX3
  for (size_t offset = 0;; offset += incr) {
    auto in{LoadU(kFullWidth8, alignment + offset)};                                     // load up to 64 bytes
    StoreU(LoadU(kFullWidth8, bases_in + offset), kFullWidth8, consensus_out + offset);  // Copy the bases
    auto mask{Ne(in, Zero(kFullWidth8))};  // Compare the input with the base, 0 means a match
    if (AllFalse(kFullWidth8, mask)) {     // true if all bases match
      nr_bases += incr;
    } else {
      nr_bases += FindKnownFirstTrue(kFullWidth8, mask);  // Add index of first non-match
      break;                                              // break out of the for loop, we're done
    }
  }
  return static_cast<int>(nr_bases);
}

// Finds the leading run of zero-valued bytes (edlib match encoding) using SIMD comparison,
// then broadcasts the specified quality character to matching positions in the output buffer.
// Uses trailing zero count on the inverted match mask to determine the run length efficiently.
HWY_ATTR int CopyMatchingQualSIMD(const uint8_t* alignment, uint8_t* quality_out, char base) {
  size_t nr_bases{0};
  const auto match_quality = Set(kFullWidth8, base);
  const auto incr{Lanes(kFullWidth8)};  // this will be 64 in AVX3
  for (size_t offset = 0;; offset += incr) {
    auto in{LoadU(kFullWidth8, alignment + offset)};           // load up to 64 bytes
    StoreU(match_quality, kFullWidth8, quality_out + offset);  // Copy the quality
    auto mask{Ne(in, Zero(kFullWidth8))};                      // Compare the input with the base, 0 means a match
    if (AllFalse(kFullWidth8, mask)) {                         // true if all bases match
      nr_bases += incr;
    } else {
      nr_bases += FindKnownFirstTrue(kFullWidth8, mask);  // Add index of first non-match
      break;                                              // break out of the for loop, we're done
    }
  }
  return static_cast<int>(nr_bases);
}

/// To find markers in the sequence, we need to search in numerous candidate locations. I have redefined the
/// LUT search (which used to use hash tables) as two subsequent LUT operations ("cascaded LUTs"). Because the
/// new approach uses a straightforward memory lookup, we now can vectorize the search and use a gather operation
/// to read multiple values from memory in one go. This is a significant improvement over the previous approach.
/// This function performs a cascaded LUT approach to find markers in the sequence. The function will return the
/// best marker that was found.
/// Assumptions: 1) num_offsets is a multiple of 16 to allow for efficient vectorization on AVX-512 platforms.
///              2) LUTs are correctly initialized.
///              3) No error checking is performed on the input data, assuming everything is OK

HWY_ATTR bool FindMarkerSIMD(const xoos::demux::CascadedLUTs& luts, xoos::demux::FixedReadRecord& record,
                             const int64_t* offsets, const int64_t* masks, const int64_t* barcode_types) {
  auto* results = reinterpret_cast<int32_t*>(record.trim_info_duplex.unfiltered_matches.data());
  auto& nr_results = record.trim_info_duplex.unfiltered_matches_count;

  const auto lut0_mask = Set(kFullWidth64S, xoos::demux::CascadedLUTs::kLut0Mask);
  const auto presence_mask = Set(kFullWidth32S, 0x0300);  // to mask in those bits that denote edist
  const auto result_mask = Set(kFullWidth32S, 0xffff);    // to mask in the 16-bit results
  const auto all_empty = Set(kFullWidth32S, 0x0300);      // 3 is the value for an empty result
  const auto multiplier =
      Set(kFullWidth64S, xoos::demux::DuplexMatch::BarcodeType::kUnknown);  // used to multiply the indices by 4

  // We expect to analyze a total of 8 offsets; only AVX-512 can handle this many offsets in one go. For the
  // other likely target, AVX2, we need to split the processing into two parts.
  for (size_t i = 0; i < 8; i += Lanes(kFullWidth64S)) {
    // Load the data from the sequence using a gather operation. This will load up to 8 positions in one go.
    auto indices{ShiftRight<2>(LoadU(kFullWidth64S, offsets + i))};
    const auto data_masks{LoadU(kFullWidth64S, masks + i)};
    const auto types{LoadU(kFullWidth64S, barcode_types + i)};
    const auto first_mask{And(lut0_mask, data_masks)};

    // The indices register have been corrected for the 2-bit offset, but now we need to re-establish the
    // correct offset, so we need to multiply the index by 4 again (note that this implicitly forces the
    // position to be a multiple of 4).
    auto indices0{ShiftLeft<2>(indices)};
    auto indices1{Add(indices0, Set(kFullWidth64S, 1))};
    auto in_raw_64_0{GatherOffset(kFullWidth64S, reinterpret_cast<const int64_t*>(record.TwoBitsSeq()), indices)};

    for (int s = 0; s < 4; ++s) {
      // We now have up to eight 64-bit inputs. Each of these inputs will provide us with markers at 8 different
      // positions (with 32 bases per position, it is easy to obtain 8 different positions for the worst-case
      // marker length and will work fine up to a marker length of 24 bases).
      // Obtain the input for the second (out of 8 positions) by shifting the input to the right by 1 base (2 bits).
      auto in_raw_64_1{ShiftRight<2>(in_raw_64_0)};

      // We now reduce the number of bases in the input to be able to use the cascaded LUTs. After the masking
      // operation, we will switch to "32-bit mode" to be able to use the cascaded LUTs.
      auto in_0_0 = And(in_raw_64_0, first_mask);
      auto in_1_0 = And(in_raw_64_1, first_mask);

      // To avoid conflict between markers, we used different entries for each marker. Take this now into
      // account by multiplying the index by the number of markers used.
      in_0_0 = Mul(in_0_0, multiplier);
      in_1_0 = Mul(in_1_0, multiplier);

      // And add the marker type so that we have loaded the correct offset for the marker.
      in_0_0 = Add(in_0_0, types);
      in_1_0 = Add(in_1_0, types);

      // Now interleave the two lower-32 bit values of the two bit inputs. I don't know the magic Highway
      // chant, but a simple shift + OR will do the trick.
      auto in_32_0 = BitCast(kFullWidth32S, Or(in_0_0, ShiftLeft<32>(in_1_0)));

      // We'll do the same for the positions.
      auto indices_32_0 = BitCast(kFullWidth32S, Or(indices0, ShiftLeft<32>(indices1)));
      // Move the positions into the high 16 bits of the 32-bit register.
      indices_32_0 = ShiftLeft<16>(indices_32_0);

      // We now have up to 16 values as input for the cascaded LUTs. Gather the contents of the first LUT; we need
      // to use 32-bit indices and therefore we use the "index" version of gather.
      auto in_lut0{GatherIndex(kFullWidth32S, reinterpret_cast<const int32_t*>(luts.Lut0().data()), in_32_0)};

      // While that request is pending, we need to prepare for the next step, i.e. gathering the offsets into the
      // second LUT.
      auto in_0_1{
          BitCast(kFullWidth64S, ShiftRight<xoos::demux::CascadedLUTs::kLut0Bits>(And(data_masks, in_raw_64_0)))};
      auto in_1_1{
          BitCast(kFullWidth64S, ShiftRight<xoos::demux::CascadedLUTs::kLut0Bits>(And(data_masks, in_raw_64_1)))};

      auto offsets_32_1{BitCast(kFullWidth32S, Or(in_0_1, ShiftLeft<32>(in_1_1)))};

      // Calculate the offsets into the second LUTs.
      offsets_32_1 = Add(offsets_32_1, in_lut0);
      // Convert to byte offsets (we now have 16-bit offsets to index into a short array).
      offsets_32_1 = ShiftLeft<1>(offsets_32_1);

      // Now gather the contents of the second LUT. Note that AVX512 only handles gathering of 32-bit and 64-bit values,
      // so we read the LUTs in 32-bit mode. This means that of the 32-bit output, only the lower 16 bits are relevant.
      // This is great, because we then can use the upper 16-bits to store an offset where the marker was found.
      auto out_0 = GatherOffset(kFullWidth32S, reinterpret_cast<const int32_t*>(luts.Lut1().data()), offsets_32_1);
      auto final_out = Or(indices_32_0, And(out_0, result_mask));

      // Which markers were found? Create a 32-bit mask that will be used to determine the presence of a marker.
      const auto masked_out = And(final_out, presence_mask);
      auto presence{Ne(masked_out, all_empty)};
      nr_results += CompressStore(final_out, presence, kFullWidth32S, results + nr_results);

      // Prepare for next iteration.
      in_raw_64_0 = ShiftRight<4>(in_raw_64_0);         // shift input by two bases
      indices0 = Add(indices0, Set(kFullWidth64S, 2));  // add two to the positions
      indices1 = Add(indices1, Set(kFullWidth64S, 2));
    }
  }

  return nr_results > 0;
}

// This function is used to find the hairpin in the sequence. The original method was developed by Igor Mandric
// (Sample demultiplexing for Hairpin Duplex Sequencing Data, 2024) and is using two sliding windows at a fixed
// distance from each other (the distance being the loop distance, which is 7 bases); the "3p" and "5p" windows
// are compared to each other and the number of matching bases is counted. If the number of matches exceeds a
// threshold, the hairpin is considered to be found.
//
// This implementation uses a similar approach, but there are some differences:
// 1) Instead of 4-bit encoding, we use 2-bit encoding. A 64-bit register will therefore contain 32 bases (instead of
//    16), which allows us to compare more bases (23 if my reasoning is correct) and likely get a more accurate result.
//    The use of 2-bit data makes the implementation more complex as we need to process 4 positions at a time.
// 2) To make the method less sensitive to indels, we are using the "Shifted Hamming Distance" (SHD) as the metric.
//    (see
//    https://academic.oup.com/bioinformatics/article/31/10/1553/176591/Shifted-Hamming-distance-a-fast-and-accurate).
//    Calculating the SHD for all four positions involves multiple shifts and XOR operations.
// 3) The method is vectorized, which means that we can process up to 8x4 = 32 positions in one go starting at the
//    specified offset position (which is assumed to be a multiple of 4).

HWY_ATTR int64_t FindHairpinSIMD(xoos::demux::FixedReadRecord& record, int64_t offset) {
  // I tried various ways to implement this function, but concluded that readability is more important than eking
  // out the last bit of performance. The function is complex enough as it is, so I will implement it in a way that
  // should be readable and maintainable.
  // We will loop over four relative positions and determine the SHD for each of these positions in parallel for
  // eight values at at time in case of AVX512; this function therefore calculates the SHD for 32 positions in one go.

  auto step{static_cast<int64_t>(Lanes(HWY_NAMESPACE::kFullWidth64) << 2)};
  offset -= std::min(offset, step);

  // We'll write 32 results (16 bits each) into this register.
  auto match_count{Zero(kFullWidth64)}, match_count_shd{Zero(kFullWidth64)};
  constexpr int kMatchShift[] = {0, 16, 32, 48};

  for (int64_t rel = 0; rel < 4; ++rel) {
    // The input offset is the first position for which we will read the "3p" part of the hairpin.
    const int64_t offset_3p = offset + rel;
    // Set up an offsets array so we can read the input data in parallel; each 64-bit value is offset by 1 byte,
    // i.e. 4 bases.
    alignas(64) constexpr int64_t kOffsets[] = {0l, 1l, 2l, 3l, 4l, 5l, 6l, 7l};
    const uint8_t* start_3p_ptr{(offset_3p >> 2) + record.TwoBitsSeq()};
    const auto in_3p_raw{
        GatherOffset(kFullWidth64S, reinterpret_cast<const int64_t*>(start_3p_ptr), Load(kFullWidth64S, kOffsets))};
    // Now we need to shift the 3p part by up to 3 bases (6 bits) to get the 2bit encoded data in the correct position.
    const int base_shift_3p = static_cast<int>(offset_3p & 3);
    const auto in_3p = ShiftRightSame(in_3p_raw, base_shift_3p + base_shift_3p);

    // The 5p part is trickier. Assuming a loop length of 7 bases and 2-bit encoding of the bases, the relative
    // position of the 5p part is 39 (32+7) bases "earlier" compared to the 3p part. But notice the following: if the
    // 3p part is perfectly aligned (base_shift_3p is 0), the 3p part would not be (because of the 39 difference in
    // offset, which is not a multiple of 4). To make sure that we always have maximum overlap between the 3p and 5p
    // parts, we therefore need an offset of 35 (28+7) and shift "in the other direction":
    // Instead of shifting 3 bases to the right, we shift 1 base to the left.
    // Instead of shifting 2 bases to the right, we shift 2 bases to the left.
    // Instead of shifting 1 base to the right, we shift 3 bases to the left.
    // Instead of shifting 0 bases to the right, we shift 4 bases to the left.
    const int64_t offset_5p = offset_3p - 28 - 7;
    const uint8_t* start_5p_ptr{(offset_5p >> 2) + record.TwoBitsSeq()};
    const auto in_5p_raw{
        GatherOffset(kFullWidth64S, reinterpret_cast<const int64_t*>(start_5p_ptr), Load(kFullWidth64S, kOffsets))};
    const int base_shift_5p = static_cast<int>(4L - (offset_5p & 3));
    auto in_5p = ShiftLeftSame(in_5p_raw, base_shift_5p + base_shift_5p);

    // Now reverse the 5p read. Use a special AVX512 instruction to do that efficiently.
#if HWY_TARGET <= HWY_AVX3_DL
    // The fastest way to reverse the 64-bit value is to use the Galois affine transformation.
    // This is a new feature only available on AVX512 processors. I have found the magic chant that allows me to use
    // this feature in Highway and also take into account that the reverse of the binary sequence 10 10 would be 10 10,
    // not 01 01 (i.e. one needs to reverse every 2 bits, not every bit).

    in_5p = hwy::HWY_NAMESPACE::BitCast(kFullWidth64S, hwy::HWY_NAMESPACE::detail::GaloisAffine(
                                                           hwy::HWY_NAMESPACE::BitCast(kFullWidth8, in_5p),
                                                           hwy::HWY_NAMESPACE::Set(kFullWidth64, 0x4080102004080102u)));
#else
    // The Galois instructions for AVX-512 are a lifesaver, but we need to implement a reverse for non-AVX512. While
    // I initially thought that this would be hard to do, I realized I can use a combination of shifts and masks and
    // a reverse lanes to do this.
    const auto mask0{Set(kFullWidth8, 0b00000011)};  // Mask for first base
    const auto mask1{Set(kFullWidth8, 0b00001100)};  // Mask for second base
    const auto mask2{Set(kFullWidth8, 0b00110000)};  // Mask for third base
    const auto mask3{Set(kFullWidth8, 0b11000000)};  // Mask for fourth base
    const auto in_5p_8bit{hwy::HWY_NAMESPACE::BitCast(kFullWidth8, in_5p)};
    auto in_5p_masked0 = ShiftLeft<6>(And(in_5p_8bit, mask0));         // shift the first base to the highest 2 bits
    const auto in_5p_masked1 = ShiftLeft<2>(And(in_5p_8bit, mask1));   // shift the second base
    auto in_5p_masked2 = ShiftRight<2>(And(in_5p_8bit, mask2));        // shift the third base
    const auto in_5p_masked3 = ShiftRight<6>(And(in_5p_8bit, mask3));  // shift the fourth base
    in_5p_masked0 = Or(in_5p_masked0, in_5p_masked1);                  // combine the first and second base
    in_5p_masked2 = Or(in_5p_masked2, in_5p_masked3);                  // combine the third and fourth base
    in_5p_masked0 = Or(in_5p_masked0, in_5p_masked2);                  // combine all bases
    in_5p = hwy::HWY_NAMESPACE::BitCast(kFullWidth64S, in_5p_masked0);
#endif
    in_5p = ReverseLaneBytes(in_5p);  // Reverse the order of the bytes per 64-bit values
    // 3p and 5p should now be perfectly aligned, both starting at position 0. We now will calculate the SHD; we'll
    // use 3p as the reference and perform shifting of the 5p part. Calculate the XOR value for each of the bits;
    // the result will be 0b11 for every matching base. By shifting the xor result by 1 base and AND-ing the result,
    // we'll get 1 bit per base that represents a match. We'll mask out the invalid bits later when we calculate the
    // popcnt().
    auto xor_result{Xor(in_3p, in_5p)};
    xor_result = And(xor_result, ShiftRight<1>(xor_result));  // this is the bit mask for rel offset 0
    auto xor_result_shd = xor_result;

    // Now calculate the SHD for other relative positions of the 5p part.
    auto shift_1p = Xor(in_3p, ShiftRight<2>(in_5p));
    shift_1p = And(shift_1p, ShiftRight<1>(shift_1p));  // this is the bit mask for rel offset 1
    xor_result_shd = Or(xor_result_shd, shift_1p);

    auto shift_1m = Xor(in_3p, ShiftLeft<2>(in_5p));
    shift_1m = And(shift_1m, ShiftRight<1>(shift_1m));  // this is the bit mask for rel offset -1
    xor_result_shd = Or(xor_result_shd, shift_1m);
    // We now have the SHD for all relative positions and need to calculate the popcnt. Most of the 3p data is valid
    // (worst-case bases 29-31 are invalid), but 5p has more invalid data because of the various shifts. Assume
    // that the first 2 bases are invalid (because of the left shifts) and that the last 7 bases are invalid. We
    // therefore can use bases 2-25 for the popcnt calculation.
    constexpr uint64_t kPopcntMask = 0x0001555555555550L;  // 23 bits are valid, max popcnt is 23.
    auto popcnt_result =
        PopulationCount(hwy::HWY_NAMESPACE::BitCast(kFullWidth64, And(xor_result, Set(kFullWidth64S, kPopcntMask))));
    auto popcnt_result_shd = PopulationCount(
        hwy::HWY_NAMESPACE::BitCast(kFullWidth64, And(xor_result_shd, Set(kFullWidth64S, kPopcntMask))));

    match_count = Or(match_count, ShiftLeftSame(popcnt_result, kMatchShift[rel]));
    match_count_shd = Or(match_count_shd, ShiftLeftSame(popcnt_result_shd, kMatchShift[rel]));
  }

  auto match_count_16{hwy::HWY_NAMESPACE::BitCast(kFullWidth16, match_count)};
  auto match_count_16_shd{hwy::HWY_NAMESPACE::BitCast(kFullWidth16, match_count_shd)};
  match_count_16 = Or(match_count_16, ShiftLeft<8>(match_count_16_shd));
  StoreU(match_count_16, kFullWidth16, reinterpret_cast<uint16_t*>(record.match_values) + offset);

  return offset;
}

}  // namespace HWY_NAMESPACE
#if HWY_ONCE

// This macro declares a static array used for dynamic dispatch.
HWY_EXPORT(To2BitSIMD);
HWY_EXPORT(FindSpaceSIMD);
HWY_EXPORT(FindNonTrivialSpaceSIMD);
HWY_EXPORT(TransformSequenceSIMD);
HWY_EXPORT(BuildPeqSIMD);
HWY_EXPORT(FindSymmetryPositionSIMD);
HWY_EXPORT(FindMarkerSIMD);
HWY_EXPORT(FindHairpinSIMD);
HWY_EXPORT(ReverseCopySIMD);
HWY_EXPORT(CopyMatchingSeqSIMD);
HWY_EXPORT(CopyMatchingQualSIMD);

void simd::TransformSequence(const char* const buf, int begin, int end, unsigned char* digit, bool reverse,
                             bool invert) {
  return HWY_DYNAMIC_DISPATCH(TransformSequenceSIMD)(reinterpret_cast<const unsigned char*>(buf), begin, end, digit,
                                                     reverse, invert);
}

void simd::ReverseCopy(const unsigned char* p_start, const unsigned char* p_end, unsigned char* p_out) {
  return HWY_DYNAMIC_DISPATCH(ReverseCopySIMD)(p_start, p_end, p_out);
}

void simd::BuildPeq(const unsigned char* digit, int length, uint64_t* bit) {
  return HWY_DYNAMIC_DISPATCH(BuildPeqSIMD)(digit, length, bit);
}

void simd::ConvertTo2Bit(const uint8_t* const p_seq, uint32_t start, uint32_t length, uint8_t* p_output) {
  return HWY_DYNAMIC_DISPATCH(To2BitSIMD)(p_seq + start, length, p_output);
}

uint64_t simd::FindFirstSpace(const uint8_t* const buf, const uint64_t begin, const uint64_t end) {
  return HWY_DYNAMIC_DISPATCH(FindSpaceSIMD)(buf, begin, end);
}

uint64_t simd::FindFirstNonTrivialSpace(const uint8_t* const buf, const uint64_t begin, const uint64_t end) {
  return HWY_DYNAMIC_DISPATCH(FindNonTrivialSpaceSIMD)(buf, begin, end);
}

int simd::CopyMatchingSeq(const unsigned char* alignment, char* bases_in, char* consensus_out) {
  return HWY_DYNAMIC_DISPATCH(CopyMatchingSeqSIMD)(alignment, reinterpret_cast<uint8_t*>(bases_in),
                                                   reinterpret_cast<uint8_t*>(consensus_out));
}

int simd::CopyMatchingQual(const unsigned char* alignment, char* quality_out, char base) {
  return HWY_DYNAMIC_DISPATCH(CopyMatchingQualSIMD)(alignment, reinterpret_cast<uint8_t*>(quality_out), base);
}

namespace xoos::demux {
int FindSymmetryPosition(FixedReadRecord& record) { return HWY_DYNAMIC_DISPATCH(FindSymmetryPositionSIMD)(record); }

bool FindMarker(const CascadedLUTs& luts, FixedReadRecord& record, const int64_t offsets[8], const int64_t* masks,
                const int64_t* types) {
  return HWY_DYNAMIC_DISPATCH(FindMarkerSIMD)(luts, record, offsets, masks, types);
}

int64_t FindHairpinSliding(FixedReadRecord& record, int64_t offset) {
  return HWY_DYNAMIC_DISPATCH(FindHairpinSIMD)(record, offset);
}

}  // namespace xoos::demux

#endif  // HWY_ONCE
