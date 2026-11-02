#pragma once

#include <iterator>
#include <memory>
#include <mutex>
#include <string>

#include <htslib/faidx.h>

#include <xoos/types/int.h>

namespace xoos::io {

struct FaidxDeleter {
  void operator()(faidx_t* const faidx) const {
    fai_destroy(faidx);
  }
};

using FaidxPtr = std::unique_ptr<faidx_t, FaidxDeleter>;

class FastaReader {
 public:
  // Holds a std::mutex, so the reader is neither copyable nor movable.
  explicit FastaReader(const std::string& fasta_file_path);
  FastaReader(const FastaReader&) = delete;
  FastaReader(FastaReader&&) = delete;
  ~FastaReader();
  FastaReader& operator=(const FastaReader&) = delete;
  FastaReader& operator=(FastaReader&&) = delete;

  /// Fetch the reference sequence for [start, end) using 0-based half-open coordinates.
  std::string GetSequence(const std::string& chrom, s64 start, s64 end);
  std::string GetSequence(const std::string& chrom, s64 start, s64 end, bool uppercase);
  /// Fetch the entire reference sequence for a chromosome.
  std::string GetSequence(const std::string& chrom);
  std::string GetSequence(const std::string& chrom, bool uppercase);
  s64 GetContigLength(const std::string& chrom);
  bool HasContig(const std::string& chrom);

  class Iterator {
   public:
    using value_type = std::pair<std::string, s32>;     // NOLINT
    using reference = value_type&;                      // NOLINT
    using pointer = value_type*;                        // NOLINT
    using iterator_category = std::input_iterator_tag;  // NOLINT
    using difference_type = std::ptrdiff_t;             // NOLINT

    Iterator(const FastaReader* reader, s32 index);

    value_type operator*() const;

    Iterator& operator++();

    bool operator!=(const Iterator& other) const;

   private:
    const FastaReader* _reader;
    s32 _index;
  };

  Iterator begin() const;  // NOLINT
  Iterator end() const;    // NOLINT

 private:
  FaidxPtr _faidx;
  std::mutex _mutex;
};

}  // namespace xoos::io
