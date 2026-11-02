#pragma once

#include <array>
#include <vector>

#include <xoos/types/int.h>

namespace xoos::demux {

class SeqMatcher;
struct DuplexMetrics;

constexpr u8 kLast2BitsMask{0x3};
constexpr u8 k00{0};
constexpr u8 k01{1};
constexpr u8 k10{2};
constexpr u8 kRightShift2{2};
// Each DNA base is encoded in 2 bits, so a base-count is converted to a bit-width by multiplying by this.
constexpr size_t kBitsPerBase{2};

/**
 * A bundle of LUTs to be used in trimming duplex data. These LUTs are meant
 * to be loaded from disk, but they might be created in memory for tests cases.
 *
 * For fastest demux performance, we actually use a cascade of LUTs; the classes used for this
 * functionality are also declared in this header file.
 */

// To minimize memory usage, we'd like to limit the amount of memory used by this structure to 16 bits.
struct DuplexMatch {
  static constexpr u32 kUninitialized{3};

  // Fixed underlying type (u16) matches the barcode_type bit-field storage: every 3-bit value is in range,
  // so casting the raw bit-field to this enum is well-defined and comparisons stay enum-to-enum (MISRA M23_061).
  // Kept as an unscoped enum (not enum class): the enumerators index fixed-size std::arrays and drive loop
  // bounds across the duplex adapters, which relies on implicit conversion to the underlying integer.
  enum BarcodeType : u16 {  // NOSONAR(S3642)
    kSID5p,
    kSID3p,
    kUMI5p,
    kUMI3p,
    kUnknown
  };

  // Length encodings for the code match. They are defined such that the value of the enum can be used to determine
  // to update the LUT table if needed.
  // longer matches are better, assign smallest value to longest match
  static constexpr u16 kPlus2{0b000};
  static constexpr u16 kPlus1{0b001};
  static constexpr u16 kZero{0b010};
  static constexpr u16 kMinus1{0b011};
  static constexpr u16 kMinus2{0b101};
  static constexpr u16 kNone{0b111};

  static constexpr auto kPos2{2};
  static constexpr auto kPos1{1};
  static constexpr auto kNeg1{-1};
  static constexpr auto kNeg2{-2};

  /**
   * @brief Maps a signed base-count delta to the corresponding Length encoding.
   *
   * @param delta The difference between the observed and nominal barcode length (-2..+2).
   * @return The encoded length value (kPlus2..kMinus2), or kZero for out-of-range deltas.
   */
  static constexpr u16 LengthFromDelta(const s32 delta) {
    switch (delta) {
      case kPos2:
        return kPlus2;
      case kPos1:
        return kPlus1;
      case 0:
        return kZero;
      case kNeg1:
        return kMinus1;
      case kNeg2:
        return kMinus2;
      default:
        return kZero;
    }
  }

  bool operator<(const DuplexMatch& other) const {
    return edist < other.edist || (edist == other.edist && length < other.length);
  }

  DuplexMatch() = default;

  u16 barcode_id : 8 {0};  // NOSONAR cpp:S2806 - u16 packs DuplexMatch into 16 bits; max 256 IDs
  u16 edist : 2 {3};       // NOSONAR cpp:S2806 - u16 packs DuplexMatch into 16 bits; edit distance 0-3
  u16 length : 3 {kNone};  // NOSONAR cpp:S2806 - u16 packs DuplexMatch into 16 bits; see Length enum
  // future extension for umi3p, umi5p, sid3p, sid5p, ... (up to 8 types)
  u16 barcode_type : 3 {kUnknown};  // NOSONAR cpp:S2806 - u16 packs DuplexMatch into 16 bits
};

// Helper struct used to store match results.
struct DuplexResult {
  DuplexResult& operator=(const DuplexResult& match) = default;

  DuplexMatch match;
  u16 pos{0};

  bool operator<(const DuplexResult& other) const {
    return match < other.match;
  }
};

// size_t min_index = std::min({ min_start, min_sid_5p, min_sid_3p, min_stop });
// I have examined how total memory usage varies with edit distance and size of the first LUT. It turns out that
// a size of 11-12 is clearly the best choice for memory usage and performance.
class CascadedLUTs {
 public:
  static constexpr size_t kNrBasesLUT0{10};  // Store 11 bases in the first LUT

  CascadedLUTs(const SeqMatcher& seqs_5p,
               const SeqMatcher& seqs_3p,
               DuplexMatch::BarcodeType barcode_type_5p,
               DuplexMatch::BarcodeType barcode_type_3p);

  static constexpr size_t kLut0Bits{kNrBasesLUT0 +
                                    kNrBasesLUT0};  // hard-coded, see table above - 9 bases, 2 bits per base
  static constexpr size_t kLut0Mask{(1 << kLut0Bits) - 1};

  // Return an approximation of the amount of memory used by the LUTs.
  size_t MemoryUsage() const {
    return (_lut0.size() * sizeof(u32)) + (_lut1.size() * sizeof(DuplexMatch)) + sizeof(CascadedLUTs);
  }

  const std::vector<u32>& Lut0() const {
    return _lut0;
  }

  const std::vector<DuplexMatch>& Lut1() const {
    return _lut1;
  }

  size_t MaxLength(const DuplexMatch::BarcodeType type) const {
    return _max_length[type];
  }

  // Given a barcode type and a match, return the length of the barcode.
  size_t Length(DuplexMatch::BarcodeType type, const DuplexMatch& match) const {
    switch (match.length) {
      case DuplexMatch::kPlus2:
        return _nominal_length[type] + 2;
      case DuplexMatch::kPlus1:
        return _nominal_length[type] + 1;
      case DuplexMatch::kZero:
        return _nominal_length[type];
      case DuplexMatch::kMinus1:
        return _nominal_length[type] - 1;
      case DuplexMatch::kMinus2:
        return _nominal_length[type] - 2;
      default:
        return 0;
    }
  }

 private:
  void AddToLUTs(const SeqMatcher& seqs, DuplexMatch::BarcodeType barcode_type);
  void AddMemory(size_t index, size_t nr_entries);

  // Inserts a single (possibly wildcard-expanded) hash entry into the cascaded LUTs, allocating LUT1 memory
  // on demand and keeping the best match per slot. Throws if two different SIDs collide on the same slot.
  void InsertHashEntry(size_t hash_value_wildcard,
                       size_t original_hash_value,
                       size_t length,
                       size_t lut1_size,
                       const DuplexMatch& lut_entry,
                       u16 barcode_id,
                       DuplexMatch::BarcodeType barcode_type);
  std::vector<u32> _lut0;                                     // first LUT, contains offset into _lut2
  std::vector<DuplexMatch> _lut1;                             // second LUT, contains the matches
  std::array<size_t, DuplexMatch::kUnknown> _nominal_length;  // nominal length of the barcodes
  std::array<size_t, DuplexMatch::kUnknown> _max_length;      // maximum length of the barcodes
};

}  // namespace xoos::demux
