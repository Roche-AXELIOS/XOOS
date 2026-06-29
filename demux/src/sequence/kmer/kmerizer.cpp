#include "kmerizer.h"

#include <xoos/util/sequence-functions.h>

#include <limits>

namespace xoos::demux::kmer {

Kmerizer::Kmerizer(std::string_view seq, size_t k)
    : _seq(seq), _k(k), _mask((1ULL << k * 2) - 1), _shift((k - 1) * 2) {}

Kmerizer::Iterator::Iterator(const Kmerizer* kmerizer, size_t pos)
    : _kmerizer(kmerizer), _pos(pos), _sub_str_len(0), _fw_kmer(0), _rv_kmer(0) {
  Next();
}

// raw hash k-mers to be compared and hashed
Kmerizer::Iterator::value_type Kmerizer::Iterator::operator*() const { return PairedKmer(_fw_kmer, _rv_kmer); }

bool Kmerizer::Iterator::operator==(const Iterator& other) const { return _pos == other._pos; }

bool Kmerizer::Iterator::operator!=(const Iterator& other) const { return !(*this == other); }

Kmerizer::Iterator& Kmerizer::Iterator::operator++() {
  Next();
  return *this;
}

void Kmerizer::Iterator::Next() {
  Step();
  while (_sub_str_len < _kmerizer->_k && _pos != std::numeric_limits<u64>::max()) {
    Step();
  }
}

void Kmerizer::Iterator::Step() {
  if (_pos < _kmerizer->_seq.size()) {
    const u8 c = sequence::kBaseToBin[static_cast<u8>(_kmerizer->_seq.at(_pos++))];
    if (c < 4) {
      _fw_kmer = (_fw_kmer << 2 | c) & _kmerizer->_mask;
      _rv_kmer = _rv_kmer >> 2 | (static_cast<BinaryKmer>(3 - c) << _kmerizer->_shift);
      if (++_sub_str_len >= _kmerizer->_k) {
        return;
      }
    } else {
      _sub_str_len = 0;
      _fw_kmer = 0;
      _rv_kmer = 0;
    }
  } else {
    _pos = std::numeric_limits<u64>::max();
  }
}

Kmerizer::Iterator Kmerizer::begin() const { return {this, 0}; }

Kmerizer::Iterator Kmerizer::end() const { return {this, std::numeric_limits<u64>::max()}; }

}  // namespace xoos::demux::kmer
