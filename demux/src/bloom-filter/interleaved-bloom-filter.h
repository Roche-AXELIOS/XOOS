#pragma once

#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

namespace xoos::demux::strand {
// An interleaved Bloom Filter is 2 bloom filters in the same bit array
// The filters are 2 equally sized bloom filters and every hash location is split into 2-bits
// Benefits: Shared header and file I/O and better cache locality if looking up both filters at the same time
//
// A Bloom filter is a probabilistic data structure for set-membership test, which has false
// positives but no false negatives. The false positive rate depends on the Bloom
// filter size and the number of hash values per set member. A larger size and,
// depending on the occupancy, more hash values, will lower the false positive rate.
//
// We enforce a power of 2 for this size to make the bit manipulation faster (modulo).

class InterleavedBloomFilter {
 public:
  InterleavedBloomFilter(size_t elements1, size_t elements2, u32 num_hashes, u32 kmer_size, f64 fpr);

  // Load an already existing interleaved bloom filter from a file
  explicit InterleavedBloomFilter(const std::filesystem::path& filename);

  void Save(const std::filesystem::path& filename) const;

  void Insert(u64 kmer, bool strand);

  // Can return false positives, but not false negatives
  std::pair<bool, bool> Contains(u64 kmer) const;

  // Total number of bits of entire combined bloom filter
  size_t Size() const;

  u32 KmerSize() const;

  std::pair<f64, f64> FalsePositiveRates() const;

  static u32 CalcOptimalNumHashes(f64 fpr, f64 occupancy);
  static u32 CalcOptimalNumHashes(f64 fpr);

 private:
  struct FileHeader {
    u32 num_hashes{0};
    u32 kmer_size{0};
    u64 size_in_bits{0};
  };

  static FileHeader LoadHeader(const std::filesystem::path& filename) {
    std::ifstream ifs(filename, std::ios::binary);
    FileHeader header;
    if (!ifs) {
      throw std::runtime_error("Could not open file: " + filename.string());
    }
    // NOSONAR(S5356,S6022) -- reinterpret_cast to char* is required to read the raw header bytes from disk;
    // binary deserialization of a trivially-copyable struct has no type-safe alternative.
    ifs.read(reinterpret_cast<char*>(&header), sizeof(header));  // NOSONAR(S5356,S6022)
    return header;
  }

  static size_t CalcArraySize(size_t num_elements, u32 num_hashes, f64 fpr);

  const FileHeader _header;
  const u64 _filter_modulo_size{(_header.size_in_bits / 2) - 1};
  std::vector<u64> _bit_array{std::vector<u64>(_header.size_in_bits / 64)};
};
}  // namespace xoos::demux::strand
