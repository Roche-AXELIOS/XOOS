#include "interleaved-bloom-filter.h"

#include <cmath>

#include "xxhash.h"

namespace xoos::demux::strand {

/**
 * @brief Create an empty Interleaved Bloom Filter (2 bloom filters interleaved together).
 * Automatically calculates the optimal size of the bloom filter (based on the maximum size of both sub-filters).
 * The interleaved bloom filter requires fewer cache misses that using two if looking up element.
 * It also allows for simpler management of the two filters (one file vs two & shared header) when serializing to disk.
 *
 * @param elements1 Number of unique elements you will be inserting for Bloom filter A.
 * @param elements2 Number of unique elements you will be inserting for Bloom filter B.
 * @param num_hashes Number of hash functions to use.
 * @param kmer_size Size of the kmer, stored for convenience but not used directly in the class.
 * @param fpr Desired false positive rate.
 */
InterleavedBloomFilter::InterleavedBloomFilter(size_t elements1, size_t elements2, u32 num_hashes, u32 kmer_size,
                                               double fpr)
    : _header(FileHeader{num_hashes, kmer_size, CalcArraySize(std::max(elements1, elements2), num_hashes, fpr) * 2}) {}

/**
 * @brief Load an Interleaved Bloom Filter from a file.
 *
 * @param filename The path to the file containing the Interleaved Bloom filter.
 * @details The file must be in the format created by the Save() method.
 */
InterleavedBloomFilter::InterleavedBloomFilter(const std::filesystem::path& filename) : _header(LoadHeader(filename)) {
  std::ifstream ifs(filename, std::ios::binary);
  ifs.seekg(0, std::ios::end);
  auto size = static_cast<std::streamoff>(ifs.tellg()) - static_cast<std::streamoff>(sizeof(_header));
  if (static_cast<u64>(size) != _header.size_in_bits / 8) {
    throw std::runtime_error("File size does not match expected size");
  }
  ifs.seekg(sizeof(_header), std::ios::beg);
  ifs.read(reinterpret_cast<char*>(_bit_array.data()), size);
}

/**
 * @brief Save the Interleaved Bloom filter to a file.
 *
 * @param filename The path to the file to save the Interleaved Bloom filter.
 */
void InterleavedBloomFilter::Save(const std::filesystem::path& filename) const {
  std::ofstream ofs(filename, std::ios::binary);
  if (!ofs) {
    throw std::runtime_error("Could not open file: " + filename.string());
  }
  ofs.write(reinterpret_cast<const char*>(&_header), sizeof(_header));

  auto size = static_cast<std::streamoff>(_bit_array.size() * sizeof(u64));
  ofs.write(reinterpret_cast<const char*>(_bit_array.data()), size);
  if (!ofs) {
    throw std::runtime_error("Could not write to file: " + filename.string());
  }
}

/**
 * @brief Inserts an item into the specified interleaved Bloom filter (Thread-Safe).
 *
 * @param kmer The kmer (64-bit number) to insert.
 * @param strand `false` targets Filter A (even bits), `true` targets Filter B (odd bits).
 */
void InterleavedBloomFilter::Insert(u64 kmer, bool strand) {
  for (uint32_t i = 0; i < _header.num_hashes; ++i) {
    // Calculate logical index [0, m-1]
    u64 logical_index = XXH64(&kmer, sizeof(kmer), static_cast<u64>(i)) & _filter_modulo_size;

    // Calculate physical index [0, 2m-1], using optimized bitwise operation equivalents
    // physical_index = 2 * logical_index + (strand ? 1 : 0)
    u64 physical_index = (logical_index << 1) | static_cast<u64>(strand);

    // Calculate bucket index and bit position using bitwise operations
    u64 bucket_index = physical_index >> 6;   // Faster than / 64
    u64 bit_in_bucket = physical_index & 63;  // Faster than % 64

    // Create the bit mask
    u64 bit_pos = 1ULL << bit_in_bucket;

    __sync_or_and_fetch(&_bit_array[bucket_index], bit_pos);
  }
}

/**
 * @brief Checks if an item is in the specified interleaved Bloom filter.
 *
 * @param kmer The kmer (64-bit number) to check.
 * @return A pair containing the result of the check for Filter A and Filter B (respectively).
 */
std::pair<bool, bool> InterleavedBloomFilter::Contains(u64 kmer) const {
  std::pair<bool, bool> result{true, true};
  for (u32 i = 0; i < _header.num_hashes; ++i) {
    // Calculate logical index [0, m-1]
    u64 logical_index = XXH64(&kmer, sizeof(kmer), static_cast<u64>(i)) & _filter_modulo_size;

    // Calculate physical index [0, 2m-1], using optimized bitwise operation equivalents
    // physical_index = 2 * logical_index
    u64 physical_index = (logical_index << 1);

    // Calculate bucket index and bit position using bitwise operations
    u64 bucket_index = physical_index >> 6;   // Faster than / 64
    u64 bit_in_bucket = physical_index & 63;  // Faster than % 64

    // Create the bit mask
    u64 bit_pos_a = 1ULL << bit_in_bucket;
    u64 bit_pos_b = 1ULL << (bit_in_bucket + 1);

    // Check Filter A (even bits) if it is not already false
    if (result.first) {
      result.first = (_bit_array[bucket_index] & bit_pos_a) != 0;
    }

    // Check Filter B (odd bits) if it is not already false
    if (result.second) {
      result.second = (_bit_array[bucket_index] & bit_pos_b) != 0;
    }

    // if both filters are false, we can skip the rest of the checks
    if (!result.first && !result.second) {
      return result;
    }
  }
  return result;
}

/**
 * @brief Estimates the False Positive Rate for each of the two interleaved filters.
 *
 * @return A pair containing the estimated FPR for Filter A (even bits) and Filter B (odd bits), respectively.
 * @details Uses the number of bits set, size of the bit vector, and the number of hash functions.
 */
std::pair<double, double> InterleavedBloomFilter::FalsePositiveRates() const {
  constexpr u64 kMaskA = 0xAAAAAAAAAAAAAAAA;  // Even bits for Filter A
  constexpr u64 kMaskB = 0x5555555555555555;  // Odd bits for Filter B

  size_t num_bits_set_a = 0;
  size_t num_bits_set_b = 0;

  for (const u64 chunk : _bit_array) {
    num_bits_set_a += std::popcount(chunk & kMaskA);
    num_bits_set_b += std::popcount(chunk & kMaskB);
  }

  size_t num_total_bits = _bit_array.size() * 64;

  // Size of a single logical filter (A or B)
  size_t num_bits_per_filter = num_total_bits / 2;

  // Calculate occupancy (fraction of bits set) for each filter
  double occupancy_a = static_cast<double>(num_bits_set_a) / static_cast<double>(num_bits_per_filter);
  double occupancy_b = static_cast<double>(num_bits_set_b) / static_cast<double>(num_bits_per_filter);

  // Estimate FPR using the approximation: FPR ≈ occupancy ^ num_hashes
  double fpr_a = pow(occupancy_a, static_cast<double>(_header.num_hashes));
  double fpr_b = pow(occupancy_b, static_cast<double>(_header.num_hashes));

  return {fpr_a, fpr_b};
}

// Return size in bits
size_t InterleavedBloomFilter::Size() const { return _header.size_in_bits; }

u32 InterleavedBloomFilter::KmerSize() const { return _header.kmer_size; }

// Typically the optimal number of hash functions is calculated based on a minimum false positive rate.
// Maximizing the entropy will minimize the false positive rate for a given number of hash functions.
// To maximize the entropy of the Bloom filter the optimal occupancy is 0.5.
// Thus, typically the best false positive rates to use are thus (0.5^num_hashes), for example:
// 1 hash function: fpr = 0.5, occupancy = 0.5
// 2 hash functions: fpr = 0.25, occupancy = 0.5
// 3 hash functions: fpr = 0.125, occupancy = 0.5
u32 InterleavedBloomFilter::CalcOptimalNumHashes(double fpr, double occupancy) {
  // check for valid fpr
  if (fpr < 0.0 || fpr > 1) {
    throw std::invalid_argument("False positive rate must be >0 and <1");
  }
  // check for valid occupancy
  if (occupancy < 0.0 || occupancy > 1) {
    throw std::invalid_argument("Occupancy must be >0 and <1");
  }

  double num_hashes = log(fpr) / log(occupancy);
  // round up to nearest integer
  num_hashes = std::round(num_hashes);
  // make sure num_hashes is at least 1, otherwise return 1
  return num_hashes > 0 ? static_cast<u32>(num_hashes) : 1;
}

/**
 * @brief Calculates the size of the bit vector needed for a single Bloom filter.
 *
 * @details The size will be rounded up to the nearest power of 2 with a minimum size of 64 bits.
 *
 * @param num_elements Number of elements to insert into the filter.
 * @param num_hashes Number of hash functions to use.
 * @param fpr Desired false positive rate.
 * @return The size of the bit vector.
 */
size_t InterleavedBloomFilter::CalcArraySize(size_t num_elements, u32 num_hashes, double fpr) {
  // This uses the exact Bloom filter formula to calculate the size of the bit array
  // The formula used in the non-approximated version, with num_elements directly used in formula in an exponent
  // This means that this version will be more accurate, though it only matters for very small values of num_elements
  // This version of formula does not allow you to solve for size_in_bits / num_elements (bits per element)
  // The size_in_bits / num_elements can be useful but was not needed we use this version of the formula
  auto hashes = static_cast<double>(num_hashes);
  auto size_in_bits = 1.0 / (1.0 - pow(1.0 - pow(fpr, 1.0 / hashes), 1.0 / hashes / static_cast<double>(num_elements)));
  // round up to nearest power of 2 and ensure minimum size of 64
  return std::bit_ceil(static_cast<size_t>(size_in_bits > 64 ? size_in_bits : 64));
}
}  // namespace xoos::demux::strand
