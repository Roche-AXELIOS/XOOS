#include "xoos/io/fasta-reader.h"

#include <mutex>
#include <stdexcept>

#include <htslib/faidx.h>

#include <xoos/error/error.h>
#include <xoos/io/malloc-ptr.h>
#include <xoos/types/int.h>
#include <xoos/util/string-functions.h>

namespace xoos::io {

FastaReader::FastaReader(const std::string& fasta_file_path) {
  const std::scoped_lock lock{_mutex};
  // Set flag to 0 to avoid automatic creation of index file
  _faidx = FaidxPtr{fai_load3(fasta_file_path.c_str(), nullptr, nullptr, 0)};
  if (_faidx == nullptr) {
    throw std::runtime_error("Unable to load fasta file or index (.fai) file");
  }
}

FastaReader::~FastaReader() {
  const std::scoped_lock lock{_mutex};
  _faidx.reset();
}

/**
 * @brief Get a sequence from a FastaReader object
 * @param chrom chromosome/contig name
 * @param start  0-based-inclusive start position
 * @param end 0-based-exclusive end position (i.e. the last position of the sequence is at end-1)
 * @param uppercase set to TRUE to convert all lowercase characters in the sequence to uppercase
 * @return sequence corresponding to (1-based region string) <chrom>:<start+1>-<end>, or 0-based half-open
 * chrom[start:end]
 */
std::string FastaReader::GetSequence(const std::string& chrom, const s64 start, s64 end, const bool uppercase) {
  std::string result;
  {
    const std::scoped_lock lock{_mutex};

    s64 len = 0;
    // NOTE: faidx_fetch_seq() is 0-based-closed, i.e. it expects that the end position is exactly the last position in
    // the desired sequence
    const MallocPtr<char> seq{faidx_fetch_seq64(_faidx.get(), chrom.c_str(), start, end - 1, &len)};
    if (seq == nullptr) {
      throw error::Error("Unable to fetch sequence '{}:{}-{}'", chrom, start + 1, end);
    }
    result = seq.get();
  }
  if (uppercase) {
    string::FastUppercase(result);
  }
  return result;
}

std::string FastaReader::GetSequence(const std::string& chrom, bool uppercase) {
  std::string result;
  {
    const std::scoped_lock lock{_mutex};
    hts_pos_t len = 0;
    const MallocPtr<char> seq{fai_fetch64(_faidx.get(), chrom.c_str(), &len)};
    if (seq == nullptr) {
      throw error::Error("Unable to fetch sequence '{}'", chrom);
    }
    result = seq.get();
  }

  if (uppercase) {
    string::FastUppercase(result);
  }
  return result;
}

FastaReader::Iterator::Iterator(const FastaReader* reader, int index) : _reader(reader), _index(index) {
}

FastaReader::Iterator::value_type FastaReader::Iterator::operator*() const {
  const char* name = faidx_iseq(_reader->_faidx.get(), _index);
  const int len = faidx_seq_len(_reader->_faidx.get(), name);
  return {std::string(name), len};
}

FastaReader::Iterator& FastaReader::Iterator::operator++() {
  ++_index;
  return *this;
}

bool FastaReader::Iterator::operator!=(const Iterator& other) const {
  return _index != other._index;
}

FastaReader::Iterator FastaReader::begin() const {
  return {this, 0};
}

FastaReader::Iterator FastaReader::end() const {
  return {this, faidx_nseq(_faidx.get())};
}

}  // namespace xoos::io
