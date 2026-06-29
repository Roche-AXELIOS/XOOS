#pragma once

#include <fmt/format.h>
#include <xoos/error/error.h>

#include <filesystem>

#include <kseq++/seqio.hpp>

#include "io/read-record.h"
#include "io/sequence-reader.h"
#include "simd/simd-functions.h"

namespace xoos::demux {

namespace fs = std::filesystem;

/**
* KZ2024 In the original implementation, KSeqppSequenceReader was a simple wrapper around klibpp::SeqStreamIn. However,
* in that approach there is considerable overhead in the allocation of strings to hold the various pieces of
* information; I therefore will convert the demux code to a new class that holds data in preallocated memory buffers
* of fixed size. So fill those buffers in an efficient way, I had to add an alternate path in kseq++; I accomplished
* this using inheritance.

* A faster alternative: avoid using strings and read directly into preallocated memory. Instead of modifying
* kseqcpp, we will inherit from it and define additional functions that allow us to read directly into buffers.
 */
template <typename TFile, typename TFunc>
class KStreamToMemory : public klibpp::KStream<TFile, TFunc, klibpp::mode::In_> {
  using CharType = klibpp::KStreamBase_::char_type;
  using SizeType = klibpp::KStreamBase_::size_type;

 public:
  using CloseType = int (*)(TFile);
  using ErrorType = std::string (*)(const TFile&);

  KStreamToMemory(TFile file, TFunc func, CloseType cfunc, ErrorType efunc, const char* filename)
      : klibpp::KStream<TFile, TFunc, klibpp::mode::In_>(file, func, cfunc), _efunc(efunc), _filename(filename) {}

 protected:
  // Because we inherit from kStream, we can use the same functions and therefore don't break backward compatibility.
  // However, we can also define new functions that allow us to read directly into buffers.
  // The first function to redefine is GetUntil, which reads directly into a buffer (code is mainly cut 'n paste,
  // with some modifications to allow for the buffer and length to be passed in).
  // No data would not copied into the output buffer if that go beyond buf_size

  inline bool GetUntil(CharType delimiter, CharType* buffer, uint& length, CharType* dret,
                       const uint buf_size)  // ks_getuntil
  {
    CharType c;
    bool gotany = false;
    if (dret) {
      *dret = 0;
    }
    SizeType i = -1;
    do {
      if (!(c = this->getc())) {
        break;
      }
      --this->begin;
      if (delimiter == this->SEP_LINE) {
        // Incorporate optimization from new seqtk (see : https://github.com/lh3/seqtk/pull/123)
        // Fabian commented on this here (https://twitter.com/kloetzl/status/1661679452479266818)
        // and suggested that std::find() may be more idiomatic.  However, I'm a bit concerned
        // that may be non-trivially slower than memchr (https://gms.tf/stdfind-and-memchr-optimizations.html).
#ifdef __GCC__
        auto* sep = reinterpret_cast<CharType*>(
            __builtin_memchr(this->buf + this->begin, '\n', this->end - this->begin));  // NOLINT
#else
        auto* sep = reinterpret_cast<CharType*>(std::memchr(this->buf + this->begin, '\n', this->end - this->begin));
#endif
        i = (sep != nullptr) ? (sep - this->buf) : this->end;
      } else if (delimiter > this->SEP_MAX) {
        for (i = this->begin; i < this->end; ++i) {
          if (this->buf[i] == delimiter) {
            break;
          }
        }
      } else if (delimiter == this->SEP_SPACE) {
        i = static_cast<SizeType>(simd::FindFirstSpace(reinterpret_cast<uint8_t*>(this->buf),
                                                       static_cast<uint64_t>(this->begin),
                                                       static_cast<uint64_t>(this->end)));
        //
        // for (i = this->begin; i < this->end; ++i) {
        //  if (IsSpace(this->buf[i])) break;
        //}
      } else if (delimiter == this->SEP_TAB) {
        i = static_cast<SizeType>(simd::FindFirstNonTrivialSpace(reinterpret_cast<uint8_t*>(this->buf),
                                                                 static_cast<uint64_t>(this->begin),
                                                                 static_cast<uint64_t>(this->end)));

        // Note: may bring this back based on future performance assessments
        // for (i = this->begin; i < this->end; ++i) {
        //   if (IsSpace(this->buf[i]) && this->buf[i] != ' ') {
        //     break;
        //   }
        // }
      } else {
        throw error::Error("Delimiter not supported {} (file: {})", delimiter, _filename);
      }

      gotany = true;
      const auto p_buf{this->buf + this->begin};
      const auto nr_items{i - this->begin};

      if (length + nr_items <= buf_size) {
        memcpy(buffer + length, p_buf, nr_items);
      }
      length += nr_items;
      this->begin = i + 1;
    } while (i >= this->end);

    if (this->err() || (this->eof() && !gotany)) {
      return false;
    }

    assert(i != -1);
    if (!this->eof() && dret) {
      *dret = this->buf[i];
    }

    // KZ024 - this test is computationally very expensive and does not seem to be necessary for Roche input files.
    return true;
  }

  void ReadSingleFixed(FixedReadRecord& rec) {
    CharType c;
    rec.Clear();  // reset all members

    this->last = false;
    if (!this->is_ready) {  // then jump to the next header line
      c = this->getc();
      if (this->fail()) {  // could be error or EOF. throw in case of error
        if (this->err()) {
          throw error::Error("Reading record start failed ({})", _efunc(this->f));
        }
        return;
      }
      if (c != '@') {
        throw error::Error("Expected '@' but got '{}' (file: {})", c, _filename);
      }
      this->is_ready = true;
    }  // else: the first header char has been read in the previous call

    uint buf_size_left = kBufferSize;
    // Read the name (description) of the sequence (usually contains instrument-related info)
    if (!this->GetUntil(this->SEP_SPACE, rec.Name(), rec.comment_offset, &c, buf_size_left)) {
      throw error::Error("Failed to read name ({})",
                         this->err() ? _efunc(this->f) : _filename + ": premature end of file");
    }
    buf_size_left = buf_size_left > rec.comment_offset ? buf_size_left - rec.comment_offset : 0;
    uint length{0};
    if (c != '\n') {  // read FASTA/Q comment
      if (!this->GetUntil(this->SEP_LINE, rec.Comment(), length, nullptr, buf_size_left)) {
        throw error::Error("Failed to read comment ({})",
                           this->err() ? _efunc(this->f) : _filename + ": premature end of file");
      }
    }
    buf_size_left = buf_size_left > length ? buf_size_left - length : 0;
    rec.seq_offset = length + rec.comment_offset;
    auto seq_offset = rec.seq_offset;
    length = 0;
    while ((c = this->getc()) && c != '+') {
      if (c == '\n') {
        continue;  // skip empty lines
      }
      rec.buf[seq_offset + length] = c;
      ++length;
      this->GetUntil(this->SEP_LINE, rec.Seq(), length, nullptr, buf_size_left);  // read the rest of the line
    }

    if (this->fail()) {
      throw error::Error("Error reading sequence bases ({})",
                         this->err() ? _efunc(this->f) : _filename + ": premature end of file");
    }

    buf_size_left = buf_size_left > length ? buf_size_left - length : 0;
    rec.qual_offset = rec.seq_offset + length;
    while ((c = this->getc()) && c != '\n') {
    }  // skip the rest of '+' line
    if (this->fail()) {  // error: no quality string
      this->is_tqs = true;
      throw error::Error("Error reading beyond sequence bases ({})",
                         this->err() ? _efunc(this->f) : _filename + ": premature end of file");
    }
    length = 0;

    while (this->GetUntil(this->SEP_LINE, rec.Qual(), length, nullptr, buf_size_left) && length < rec.SeqLen()) {
    }
    rec.end_offset = rec.qual_offset + length;
    if (this->err()) {
      this->is_tqs = true;
      throw error::Error("Error reading sequence qualities ({})", _efunc(this->f));
    }
    this->is_ready = false;               // we have not come to the next header line
    if (rec.SeqLen() != rec.QualLen()) {  // error: qual string is of a different length
      this->is_tqs = true;                // should return here
      throw error::Error("Length mismatch between sequence and quality strings (file: {})", _filename);
    }
    rec.SetStatus((rec.end_offset <= kBufferSize) ? FixedReadRecord::Status::kRead
                                                  : FixedReadRecord::Status::kTooLongFail);
    this->last = true;
    ++this->counter;
  }

 public:
  inline BatchStatistics ReadBatchFixed(FixedReadRecordBatch& batch, size_t max_batch_size) {
    batch.start_time = std::chrono::high_resolution_clock::now();
    auto start_pos{this->begin};

    batch.num_bytes = 0ul;
    batch.num_bases = 0ul;
    for (batch.num_records = 0ul; batch.num_records < max_batch_size; ++batch.num_records) {
      auto& record = (*batch.records)[batch.num_records];
      ReadSingleFixed(record);  // read fastq record
      // Now convert the sequence into 2-bit code
      record.InitTwoBitSeq();
      batch.num_bases += record.SeqLen();
      if (this->eof()) {
        batch.num_bytes += this->begin - start_pos;
        batch.end_time = std::chrono::high_resolution_clock::now();
        return {batch.num_records, batch.num_bases};
      }
    }
    batch.num_bytes += this->begin - start_pos;
    batch.end_time = std::chrono::high_resolution_clock::now();
    return {batch.num_records, batch.num_bases};
  }

 private:
  ErrorType _efunc;
  std::string _filename;
};

class SeqStreamMemoryIn : public KStreamToMemory<gzFile, int (*)(gzFile_s*, void*, unsigned int)> {
 public:
  /* Typedefs */
  typedef KStreamToMemory<gzFile, int (*)(gzFile_s*, void*, unsigned int)> BaseType;

  /* Lifecycle */
  explicit SeqStreamMemoryIn(const char* filename)
      : BaseType(
            gzopen(filename, "r"),
            // overriding gzread function with a wrapper as gzread's return code doesn't match with the base class's
            // (klibpp::KStream) expectation; when gzread returns a <= 0 value it could mean EOF or I/O errors. The
            // wrapper changes that to return 0 for EOF and -1 in case of errors.
            [](gzFile_s* file, void* buf, unsigned int len) -> int {
              auto ret = gzread(file, buf, len);
              if (ret > 0) {  // > 0: num-bytes read
                return ret;
              } else {  // <= 0: eof or error
                int err_num;
                gzerror(file, &err_num);
                return err_num == Z_OK ? 0 : -1;  // Z_OK: eof
              }
            },
            gzclose,
            // lambda to get a human-readable string describing the last error that occurred on the input gzFile and use
            // it in error reporting.
            [](const gzFile& file) -> std::string {
              int err_num;
              return gzerror(file, &err_num);
            },
            filename) {}
};

// Original KSeqppSequenceReader implementation - considerable overhead because of string copying. Will be retired.
class KSeqppSequenceReader : public SequenceReader {
 public:
  explicit KSeqppSequenceReader(const fs::path& path);

  ~KSeqppSequenceReader() override = default;

  /// Reads the next batch of records into the provided arena. It is assumed that the batch memory
  /// already has been pre-allocated and the desired number of records (i.e. batch size) is determined
  /// by the number of pre-allocated records in batch. This function updates the contents of batch and
  /// returns the number of records read (which is usually equal to the number of pre-allocated records).
  BatchStatistics ReadBatchIntoArena(FixedReadRecordBatch& batch) override;

 private:
  SeqStreamMemoryIn _seq_stream_in;
};
}  // namespace xoos::demux
