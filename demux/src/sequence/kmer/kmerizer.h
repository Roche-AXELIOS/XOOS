#pragma once

#include <xoos/types/int.h>

#include <climits>
#include <string_view>

namespace xoos::demux::kmer {
// 64-bit integer to store binary encoded k-mer, but setting it here we can consider changing it
using BinaryKmer = u64;

struct PairedKmer {
  BinaryKmer fw;
  BinaryKmer rv;
};

/**
 * @brief Class to iterate over k-mers in a sequence, and k-mer are represented in binary format.
 * The maximum k supported is 32, as each base is represented by 2 bits.
 *
 * Non-ACGT bases are treated as k-mer separators, so if a non-ACGT base is found,
 * the k-mer is reset and the next k-mer will start after the non-ACGT base.
 *
 * The iterator returns a pair of k-mers: the forward k-mer and the reverse complement k-mer.
 * To obtain the canonical k-mer do a straightforward comparison of the two k-mers numerically.
 *
 * Example:
 *   Sequence: ACGTNNACGT
 *   k = 3
 *   K-mers: ACG, CGT, ACG, CGT
 *   (Non-ACGT bases cause the search to reset)
 */
class Kmerizer {
 public:
  Kmerizer(std::string_view seq, size_t k);
  ~Kmerizer() = default;

  class Iterator {
   public:
    using value_type = PairedKmer;                      // NOLINT
    using reference = value_type&;                      // NOLINT
    using pointer = value_type*;                        // NOLINT
    using iterator_category = std::input_iterator_tag;  // NOLINT
    using difference_type = std::ptrdiff_t;             // NOLINT

    Iterator(const Kmerizer* kmerizer, size_t pos);

    value_type operator*() const;
    bool operator==(const Iterator& other) const;
    bool operator!=(const Iterator& other) const;
    Iterator& operator++();

   private:
    const Kmerizer* _kmerizer;

    // position of the current on the sequence
    // may not equal (position + 1 - k) k-mers returned because of non-ACGT bases
    size_t _pos;

    // The size of the current bases from last non-ACGT base, reset to 0 when non-ACGT base is found
    size_t _sub_str_len;

    // k-mer direction is relative to sequence
    BinaryKmer _fw_kmer;
    BinaryKmer _rv_kmer;

    void Next();
    void Step();
  };

  Iterator begin() const;  // NOLINT
  Iterator end() const;    // NOLINT

 private:
  std::string_view _seq;
  const u32 _k;
  // mask to get the last k bases
  const BinaryKmer _mask;
  // how much to shift the reverse complement bases by when appending base determined by k
  const u64 _shift;
};
}  // namespace xoos::demux::kmer
