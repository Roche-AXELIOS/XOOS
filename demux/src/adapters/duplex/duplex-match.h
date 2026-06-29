#pragma once

#include <xoos/types/int.h>

#include <array>
#include <vector>

namespace xoos::demux {

class SeqMatcher;
struct DuplexMetrics;

constexpr u8 kLast2BitsMask{0x3};
constexpr u8 k00{0};
constexpr u8 k01{1};
constexpr u8 k10{2};
constexpr u8 kRightShift2{2};

/**
 * A bundle of LUTs to be used in trimming duplex data. These LUTs are meant
 * to be loaded from disk, but they might be created in memory for tests cases.
 *
 * For fastest demux performance, we actually use a cascade of LUTs; the classes used for this
 * functionality are also declared in this header file.
 */

// KZ2024 - to minimize memory usage, we'd like to limit the amount of memory used by this structure to 16 bits.
struct DuplexMatch {
  constexpr static u32 kUninitialized{3};

  enum BarcodeType { kSID5p, kSID3p, kUMI5p, kUMI3p, kUnknown };

  // Length encodings for the code match. They are defined such that the value of the enum can be used to determine
  // to update the LUT table if needed.
  static const u16 kPlus2{0b000};  // longer matches are better, assign smallest value to longest match
  static const u16 kPlus1{0b001};
  static const u16 kZero{0b010};
  static const u16 kMinus1{0b011};
  static const u16 kMinus2{0b101};
  static const u16 kNone{0b111};

  bool operator<(const DuplexMatch& other) const {
    return edist < other.edist || (edist == other.edist && length < other.length);
  }

  DuplexMatch() = default;

  u16 barcode_id : 8 {0};           // Max 256 IDs
  u16 edist : 2 {3};                // Edit distance
  u16 length : 3 {kNone};           // See Length enum
  u16 barcode_type : 3 {kUnknown};  // future extension for umi3p, umi5p, sid3p, sid5p, ... (up to 8 types)
};

// Helper struct used to store match results.
struct DuplexResult {
  DuplexResult& operator=(const DuplexResult& match) = default;

  DuplexMatch match;
  u16 pos{0};

  bool operator<(const DuplexResult& other) const { return match < other.match; }
};

// size_t min_index = std::min({ min_start, min_sid_5p, min_sid_3p, min_stop });
// I have examined how total memory usage varies with edit distance and size of the first LUT. It turns out that
// a size of 11-12 is clearly the best choice for memory usage and performance.
class CascadedLUTs {
 public:
  static constexpr size_t kNrBasesLUT0{10};  // Store 11 bases in the first LUT

  CascadedLUTs(const SeqMatcher& seqs_5p, const SeqMatcher& seqs_3p, DuplexMatch::BarcodeType barcode_type_5p,
               DuplexMatch::BarcodeType barcode_type_3p);

  static constexpr size_t kLut0Bits{kNrBasesLUT0 +
                                    kNrBasesLUT0};  // hard-coded, see table above - 9 bases, 2 bits per base
  static constexpr size_t kLut0Mask{(1 << kLut0Bits) - 1};

  // Return an approximation of the amount of memory used by the LUTs.
  size_t MemoryUsage() const {
    return _lut0.size() * sizeof(u32) + _lut1.size() * sizeof(DuplexMatch) + sizeof(CascadedLUTs);
  }

  const std::vector<u32>& Lut0() const { return _lut0; }

  const std::vector<DuplexMatch>& Lut1() const { return _lut1; }

  size_t MaxLength(DuplexMatch::BarcodeType type) const { return _max_length[type]; }

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
  std::vector<u32> _lut0;                                     // first LUT, contains offset into _lut2
  std::vector<DuplexMatch> _lut1;                             // second LUT, contains the matches
  std::array<size_t, DuplexMatch::kUnknown> _nominal_length;  // nominal length of the barcodes
  std::array<size_t, DuplexMatch::kUnknown> _max_length;      // maximum length of the barcodes
};

}  // namespace xoos::demux
