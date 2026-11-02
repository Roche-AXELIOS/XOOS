#include "interleaved-bloom-filter.h"

#include <cmath>

#include <xoos/types/float.h>

#include "xxhash.h"

namespace xoos::demux::strand {

constexpr size_t kNumInterleavedFilters = 2;
constexpr u64 kBitsPerByte = 8;
constexpr u64 kBitsPerBucket = 64;
constexpr u64 kBucketIndexShift = 6;
constexpr u64 kBucketBitMask = kBitsPerBucket - 1;
constexpr u64 kInterleavedFilterShift = 1;
constexpr u64 kAdjacentFilterBitOffset = 1;
constexpr u64 kMinimumLogicalFilterBits = 64;
constexpr u32 kMinimumNumHashes = 1;
// The logarithmic Bloom-filter calculations require probabilities in the open interval (0, 1).
constexpr f64 kProbabilityLowerExclusiveBound = 0.0;
constexpr f64 kProbabilityUpperExclusiveBound = 1.0;
// Optimal occupancy for maximizing entropy (and thus minimizing the false positive rate).
constexpr f64 kDefaultOccupancy = 0.5;

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
InterleavedBloomFilter::InterleavedBloomFilter(
    const size_t elements1, const size_t elements2, const u32 num_hashes, const u32 kmer_size, const f64 fpr)
    : _header(FileHeader{num_hashes,
                         kmer_size,
                         CalcArraySize(std::max(elements1, elements2), num_hashes, fpr) * kNumInterleavedFilters}) {
}

/**
 * @brief Load an Interleaved Bloom Filter from a file.
 *
 * @param filename The path to the file containing the Interleaved Bloom filter.
 * @details The file must be in the format created by the Save() method.
 */
InterleavedBloomFilter::InterleavedBloomFilter(const std::filesystem::path& filename) : _header(LoadHeader(filename)) {
  std::ifstream ifs(filename, std::ios::binary);
  ifs.seekg(0, std::ios::end);
  const auto size = static_cast<std::streamoff>(ifs.tellg()) - static_cast<std::streamoff>(sizeof(_header));
  if (static_cast<u64>(size) != _header.size_in_bits / kBitsPerByte) {
    throw std::runtime_error("File size does not match expected size");
  }
  ifs.seekg(sizeof(_header), std::ios::beg);
  // NOSONAR(S5356,S6022) -- reinterpret_cast to char* is required to read the raw bit array bytes from disk;
  // binary serialization of a trivially-copyable buffer has no type-safe alternative.
  ifs.read(reinterpret_cast<char*>(_bit_array.data()), size);  // NOSONAR(S5356,S6022)
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
  // NOSONAR(S5356,S6022) -- reinterpret_cast to const char* is required to serialize the raw header bytes to disk;
  // binary serialization of a trivially-copyable struct has no type-safe alternative.
  ofs.write(reinterpret_cast<const char*>(&_header), sizeof(_header));  // NOSONAR(S5356,S6022)

  const auto size = static_cast<std::streamoff>(_bit_array.size() * sizeof(u64));
  // NOSONAR(S5356,S6022) -- reinterpret_cast to const char* is required to serialize the raw bit array bytes to disk.
  ofs.write(reinterpret_cast<const char*>(_bit_array.data()), size);  // NOSONAR(S5356,S6022)
  // NOSONAR(M23_002) -- ofs.write can fail and set the stream's failbit, so this check is runtime-dependent;
  // the analyzer does not model the stream state change from write().
  if (!ofs) {  // NOSONAR(M23_002)
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
  for (u32 i = 0; i < _header.num_hashes; ++i) {
    // Calculate logical index [0, m-1]
    // NOSONAR(S5356) -- &kmer is passed to the XXH64 C hashing API, which takes const void*; feeding the
    // raw bytes of the key to the hash function is required and has no type-safe alternative.
    const u64 logical_index = XXH64(&kmer, sizeof(kmer), static_cast<u64>(i)) & _filter_modulo_size;  // NOSONAR(S5356)

    // Calculate physical index [0, 2m-1], using optimized bitwise operation equivalents
    // physical_index = 2 * logical_index + (strand ? 1 : 0)
    const u64 physical_index = (logical_index << kInterleavedFilterShift) | (strand ? 1UL : 0UL);

    // Calculate bucket index and bit position using bitwise operations
    // (shift is faster than / 64, mask is faster than % 64)
    const u64 bucket_index = physical_index >> kBucketIndexShift;
    const u64 bit_in_bucket = physical_index & kBucketBitMask;

    // Create the bit mask
    const u64 bit_pos = 1ULL << bit_in_bucket;

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
    // NOSONAR(S5356) -- &kmer is passed to the XXH64 C hashing API, which takes const void*; feeding the
    // raw bytes of the key to the hash function is required and has no type-safe alternative.
    const u64 logical_index = XXH64(&kmer, sizeof(kmer), static_cast<u64>(i)) & _filter_modulo_size;  // NOSONAR(S5356)

    // Calculate physical index [0, 2m-1], using optimized bitwise operation equivalents
    // physical_index = 2 * logical_index
    const u64 physical_index = (logical_index << kInterleavedFilterShift);

    // Calculate bucket index and bit position using bitwise operations
    // (shift is faster than / 64, mask is faster than % 64)
    const u64 bucket_index = physical_index >> kBucketIndexShift;
    const u64 bit_in_bucket = physical_index & kBucketBitMask;

    // Create the bit mask
    const u64 bit_pos_a = 1ULL << bit_in_bucket;
    const u64 bit_pos_b = 1ULL << (bit_in_bucket + kAdjacentFilterBitOffset);

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
std::pair<f64, f64> InterleavedBloomFilter::FalsePositiveRates() const {
  // Even bits for Filter A
  constexpr u64 kMaskA = 0xAAAA'AAAA'AAAA'AAAAU;
  // Odd bits for Filter B
  constexpr u64 kMaskB = 0x5555'5555'5555'5555U;

  size_t num_bits_set_a = 0;
  size_t num_bits_set_b = 0;

  for (const u64 chunk : _bit_array) {
    num_bits_set_a += static_cast<size_t>(std::popcount(chunk & kMaskA));
    num_bits_set_b += static_cast<size_t>(std::popcount(chunk & kMaskB));
  }

  const size_t num_total_bits = _bit_array.size() * kBitsPerBucket;

  // Size of a single logical filter (A or B)
  const size_t num_bits_per_filter = num_total_bits / kNumInterleavedFilters;

  // Calculate occupancy (fraction of bits set) for each filter
  const f64 occupancy_a = static_cast<f64>(num_bits_set_a) / static_cast<f64>(num_bits_per_filter);
  const f64 occupancy_b = static_cast<f64>(num_bits_set_b) / static_cast<f64>(num_bits_per_filter);

  // Estimate FPR using the approximation: FPR ≈ occupancy ^ num_hashes
  const f64 fpr_a = std::pow(occupancy_a, static_cast<f64>(_header.num_hashes));
  const f64 fpr_b = std::pow(occupancy_b, static_cast<f64>(_header.num_hashes));

  return {fpr_a, fpr_b};
}

// Return size in bits
size_t InterleavedBloomFilter::Size() const {
  return _header.size_in_bits;
}

u32 InterleavedBloomFilter::KmerSize() const {
  return _header.kmer_size;
}

// Typically the optimal number of hash functions is calculated based on a minimum false positive rate.
// Maximizing the entropy will minimize the false positive rate for a given number of hash functions.
// To maximize the entropy of the Bloom filter the optimal occupancy is 0.5.
// Thus, typically the best false positive rates to use are thus (0.5^num_hashes), for example:
// 1 hash function: fpr = 0.5, occupancy = 0.5
// 2 hash functions: fpr = 0.25, occupancy = 0.5
// 3 hash functions: fpr = 0.125, occupancy = 0.5
u32 InterleavedBloomFilter::CalcOptimalNumHashes(const f64 fpr, const f64 occupancy) {
  // check for valid fpr
  if (std::isnan(fpr) || fpr <= kProbabilityLowerExclusiveBound || fpr >= kProbabilityUpperExclusiveBound) {
    throw std::invalid_argument("False positive rate must be >0 and <1");
  }
  // check for valid occupancy
  if (std::isnan(occupancy) || occupancy <= kProbabilityLowerExclusiveBound ||
      occupancy >= kProbabilityUpperExclusiveBound) {
    throw std::invalid_argument("Occupancy must be >0 and <1");
  }

  // round up to nearest integer
  const f64 num_hashes = std::round(std::log(fpr) / std::log(occupancy));
  // make sure num_hashes is at least 1, otherwise return 1
  return num_hashes > 0 ? static_cast<u32>(num_hashes) : kMinimumNumHashes;
}

u32 InterleavedBloomFilter::CalcOptimalNumHashes(const f64 fpr) {
  return CalcOptimalNumHashes(fpr, kDefaultOccupancy);
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
size_t InterleavedBloomFilter::CalcArraySize(const size_t num_elements, const u32 num_hashes, const f64 fpr) {
  // This uses the exact Bloom filter formula to calculate the size of the bit array
  // The formula used in the non-approximated version, with num_elements directly used in formula in an exponent
  // This means that this version will be more accurate, though it only matters for very small values of num_elements
  // This version of formula does not allow you to solve for size_in_bits / num_elements (bits per element)
  // The size_in_bits / num_elements can be useful but was not needed we use this version of the formula
  const auto hashes = static_cast<f64>(num_hashes);
  const f64 size_in_bits =
      1.0 / (1.0 - std::pow(1.0 - std::pow(fpr, 1.0 / hashes), 1.0 / hashes / static_cast<f64>(num_elements)));
  // round up to nearest power of 2 and ensure minimum size of 64
  return std::bit_ceil(
      static_cast<size_t>(size_in_bits > kMinimumLogicalFilterBits ? size_in_bits : kMinimumLogicalFilterBits));
}
}  // namespace xoos::demux::strand
