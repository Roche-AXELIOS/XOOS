#include "seq-lut.h"

#include <stdexcept>

#include "simd/simd-functions.h"

namespace xoos::demux {
SeqLut::FindResult SeqLut::SimpleFind(const std::string_view& seq) const {
  // Reroute this implementation to the new one using integer hash values - this
  // allows us to reuse current code and tests.
  return SimpleFind(seq, 0, seq.length());
}

const BarcodePool& SeqLut::Pool() const { return _pool; }

std::string SeqLut::ConvertFrom2Bit(const uint8_t* p_input, uint length) {
  constexpr char kEncoding[4]{'A', 'C', 'T', 'G'};
  const auto length2{(length + 3) / 4};
  std::string out;
  for (uint i = 0; i < length2; ++i) {
    auto c{p_input[i]};
    for (int j = 0; j < 4; ++j) {
      out.push_back(kEncoding[c & 3]);
      c >>= 2;
    }
  }

  while (out.length() > length) {
    out.pop_back();
  }
  return out;
}

uint64_t SeqLut::HashValue(const std::string_view& seq, uint start, uint length) {
  uint64_t hash{0UL};
  simd::ConvertTo2Bit(reinterpret_cast<const uint8_t*>(seq.data()), start, length, reinterpret_cast<uint8_t*>(&hash));
  return hash;
}

SeqLut::SeqLut(const IntLut& lut, const BarcodePool& pool) : _lut(lut), _pool(pool) {}  // NOLINT

SeqLut::FindResult SeqLut::SimpleFind(const std::string_view& seq, uint spos, uint epos) const {
  const auto length = epos - spos;
  return _lut(length, HashValue(seq, spos, length));
}

IntLut::IntLut(uint max_length) : _hash_tables(32) {
  if (max_length < kDefaultMaxLength) {
    throw std::invalid_argument("Must be at least 4 long");
  }
  // allocate enough memory to hold all arbitrary combinations of hash values for max length
  // Compute the prefilter mask. Assuming a maximum edit distance of 2, the prefilter length should be 2 * 2 less
  // than the maximum length (which is seq_len + max edit distance).
  // the prefilter length must be larger than 3 and in practice it will always be much larger than 3
  auto prefilter_length_2 = std::max(3u, 2 * (max_length - 4));
  // We have constraints; to optimize cache efficiency, we do need to limit the memory for the prefilter.
  // 64 KB (0xffff): 8 bases + 1 1/2 encoded in byte -> 9 1/2 bases
  // 128 KB (0x1ffff): 8 1/2 bases + 1 1/2 encoded in byte -> 10 bases
  // 512 KB (0x7ffff): 9 1/2 bases + 1 1/2 encoded in byte -> 11 bases
  // 2048 KB (0x1fffff): 10 1/2 bases + 1 1/2 encoded in byte -> 12 bases -> 24 bits
  constexpr uint kMaxNumberBits{24};
  prefilter_length_2 = std::min(prefilter_length_2, kMaxNumberBits);
  size_t nr_prefilter{1ul << (prefilter_length_2 - 3)};
  _prefilter.resize(nr_prefilter, 0);
  _prefilter_mask = static_cast<int>(nr_prefilter - 1);
}

void IntLut::Add(uint length, uint64_t hash_val, const BarcodeMatch& barcode) {
  // add to the appropriate hash-based map
  _hash_tables[length].emplace(hash_val, barcode);
  // The filter LUT supports up to 16 bases
  length = std::min(length, 16u);
  // The lowest 3 bases are stored in the LUT, one bit for every possible combination
  hash_val &= 0xffffffff;
  auto set_bits{kSetBits[hash_val & 7]};  // determine the bits that need to be set
  hash_val >>= 3;
  hash_val &= _prefilter_mask;
  _prefilter[hash_val] |= set_bits;  // the bit encodes presence in up to 11 bases
}
}  // namespace xoos::demux
