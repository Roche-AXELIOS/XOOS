#pragma once

#include <iterator>
#include <memory>
#include <mutex>
#include <string>

#include <htslib/faidx.h>

#include <xoos/types/int.h>

namespace xoos::io {

struct FaidxDeleter {
  void operator()(faidx_t* faidx) const {
    fai_destroy(faidx);
  }
};

using FaidxPtr = std::unique_ptr<faidx_t, FaidxDeleter>;

class FastaReader {
 public:
  explicit FastaReader(const std::string& fasta_file_path);
  ~FastaReader();

  std::string GetSequence(const std::string& chrom, s64 start, s64 end, bool uppercase = true);
  std::string GetSequence(const std::string& chrom, bool uppercase = true);

  class Iterator {
   public:
    using value_type = std::pair<std::string, int>;     // NOLINT
    using reference = value_type&;                      // NOLINT
    using pointer = value_type*;                        // NOLINT
    using iterator_category = std::input_iterator_tag;  // NOLINT
    using difference_type = std::ptrdiff_t;             // NOLINT

    Iterator(const FastaReader* reader, int index);

    value_type operator*() const;

    Iterator& operator++();

    bool operator!=(const Iterator& other) const;

   private:
    const FastaReader* _reader;
    int _index;
  };

  Iterator begin() const;  // NOLINT
  Iterator end() const;    // NOLINT

 private:
  FaidxPtr _faidx;
  std::mutex _mutex;
};

}  // namespace xoos::io
