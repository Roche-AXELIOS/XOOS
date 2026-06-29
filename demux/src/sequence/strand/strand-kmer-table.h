#pragma once

#include <xoos/types/int.h>

#include <atomic>
#include <vector>

namespace xoos::demux::strand {
// Strand value bits
static constexpr u64 kForwardStrandBit = 1ULL << (63);                        // 0x4000000000000000
static constexpr u64 kReverseStrandBit = 1ULL << (62);                        // 0x8000000000000000
static constexpr u64 kBothStrandBit = kForwardStrandBit | kReverseStrandBit;  // 0xC000000000000000

// Stores canonical k-mers (key) and associated flags (value) concurrently.
class StrandKmerTable {
 public:
  static constexpr u64 kEmptyEntry = 0;  // Empty entry in the table.

  static std::pair<u64, u64> DecodeEntry(u64 table_entry);

  explicit StrandKmerTable(u64 max_elements, double target_load_factor = 0.85);

  // Prevent copy and assign operations.
  StrandKmerTable(const StrandKmerTable&) = delete;
  StrandKmerTable& operator=(const StrandKmerTable&) = delete;

  void InsertOrUpdate(u64 key, bool reversed);

  const std::vector<std::atomic<u64>>& GetTable() const;

  std::size_t Size() const;

 private:
  // Using Most Significant Bits (MSB) for the 2-bit strand info value.
  static constexpr u64 kKeyMask = (1ULL << (62)) - 1;  // 0x3FFFFFFFFFFFFFFF

  // Helper to calculate table size during construction.
  static std::size_t CalculateTableSize(u64 max_elements, double target_load_factor);

  const std::size_t _table_size;
  std::vector<std::atomic<u64>> _table;
};
}  // namespace xoos::demux::strand
