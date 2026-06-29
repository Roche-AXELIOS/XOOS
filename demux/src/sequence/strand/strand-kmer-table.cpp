#include "strand-kmer-table.h"

#include <bit>        // For std::bit_ceil (C++20)
#include <cmath>      // For std::ceil
#include <limits>     // For numeric_limits checks
#include <stdexcept>  // For exceptions
#include <vector>     // Included again for implementation details if needed

#include "xxhash.h"

namespace xoos::demux::strand {

// This class is a hash table for storing canonical k-mer stand information
// Assume k-mer size is less than 31 bases, and canonical k-mer is the larger of the two complements
// The table is implemented as a vector of atomic 64-bit integers
// The key is a 62 bit value, and the value is a 2 bit strand information
// Concurrent access modification is handled using atomic and compare-and-swap (CAS) operations

/**
 * @brief Extracts the key (lower 62 bits) from a raw 64-bit table entry.
 * @param table_entry The raw 64-bit entry from the hash table.
 * @return u64 The extracted key.
 */
std::pair<u64, u64> StrandKmerTable::DecodeEntry(u64 table_entry) {
  // Extract the key (lower 62 bits) using the mask
  u64 key = table_entry & kKeyMask;

  // Extract the value (upper 2 bits) using right shift
  u64 value = table_entry & kBothStrandBit;

  return {key, value};
}

/**
 * @brief Calculate the table size as a power of 2 based on the maximum elements and target load factor.
 * @param max_elements Maximum number of elements expected in the table.
 * @param target_load_factor Desired load factor for the hash table.
 * @return The calculated table size.
 * @throws std::invalid_argument If the load factor is not between 0.0 and 1.0 or max_elements is 0.
 * @throws std::overflow_error If the calculated table size exceeds the maximum size_t value.
 */
std::size_t StrandKmerTable::CalculateTableSize(u64 max_elements, double target_load_factor) {
  if (target_load_factor <= 0.0 || target_load_factor >= 1.0) {
    throw std::invalid_argument("Load factor must be between 0.0 and 1.0");
  }
  if (max_elements == 0) {
    throw std::invalid_argument("Maximum elements must be greater than 0");
  }

  const u64 required_slots = std::ceil(static_cast<double>(max_elements) / target_load_factor);
  // Ensure required_slots is a power of 2
  const u64 table_size_p2 = std::bit_ceil(required_slots);

  if (table_size_p2 > std::numeric_limits<std::size_t>::max()) {
    throw std::overflow_error("Calculated table size exceeds size_t maximum.");
  }

  return static_cast<std::size_t>(table_size_p2);
}

/**
 * @brief Construct a StrandKmerTable with a calculated table size.
 * @param max_elements Maximum number of elements expected in the table.
 * @param target_load_factor Desired load factor for the hash table.
 */
StrandKmerTable::StrandKmerTable(u64 max_elements, double target_load_factor)
    : _table_size(CalculateTableSize(max_elements, target_load_factor)), _table(_table_size) {}

/**
 * @brief Insert or update a k-mer and its strand information in the hash table.
 * @param key 31-mer (62 bits) kmer in binary format, cannot be 0.
 * @param strand 0 for forward strand, 1 for reverse strand.
 * @throws std::runtime_error If the hash table is full and no slot is found.
 */
void StrandKmerTable::InsertOrUpdate(u64 key, bool reversed) {
  // Use bitwise AND for modulo since _table_size is a power of 2.
  const std::size_t initial_index = XXH64(&key, sizeof(key), static_cast<u64>(0)) & (_table_size - 1);
  u64 value_bits = reversed ? kForwardStrandBit : kReverseStrandBit;  // 1 forward, 0 for reverse

  for (std::size_t attempt = 0; attempt < _table_size; ++attempt) {
    // Quadratic Probing: (attempt^2 + attempt) / 2 offset.
    u64 offset = (static_cast<u64>(attempt) * attempt + attempt) / 2;
    // Use bitwise AND for power-of-2 modulo. This is faster than using % operator.
    const std::size_t current_index = (initial_index + offset) & (_table_size - 1);

    std::atomic<u64>& slot = _table[current_index];
    u64 current_entry = slot.load(std::memory_order_relaxed);

    // CAS retry loop
    while (true) {
      // Case 1: Empty slot.
      if (current_entry == 0) {
        u64 new_entry = value_bits | key;
        if (slot.compare_exchange_weak(current_entry, new_entry, std::memory_order_acq_rel,
                                       std::memory_order_acquire)) {
          return;  // Success
        }
        // CAS failed, current_entry updated, retry loop
        continue;
      }

      // Slot not empty.
      u64 entry_key = current_entry & kKeyMask;

      // Case 2: Key matches.
      if (entry_key == key) {
        slot.fetch_or(value_bits, std::memory_order_release);
        return;  // Success
      }

      // Case 3: Collision. Break inner loop to try next probe index.
      break;
    }  // End CAS retry loop
  }  // End for loop (probing attempts)

  // Should only be reached if the table is truly full and no slot was found.
  throw std::runtime_error("Hash table insert failed (table likely full).");
}

/**
 * @brief Get the raw underlying table with of atomic 64-bit integers.
 * @return A constant reference to the table.
 */
const std::vector<std::atomic<u64>>& StrandKmerTable::GetTable() const { return _table; }

/**
 * @brief Get the size of the raw hash table.
 * @return The size of the table.
 */
std::size_t StrandKmerTable::Size() const { return _table_size; }
}  // namespace xoos::demux::strand
