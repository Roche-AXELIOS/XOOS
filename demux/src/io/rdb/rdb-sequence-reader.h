#pragma once

#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

#include "io/memory-file-stream.h"
#include "io/sequence-reader.h"

namespace xoos::demux {
namespace fs = std::filesystem;

// current specification for cycle id
constexpr u32 kMinBrightCycleId = 0;
constexpr u32 kMaxBrightCycleId = 19;

/**
 * RDB types defined by run_id header
 */
enum class RdbRunIdType {
  kUnknownRdbRunIdType,
  kHtp1RdbRunIdType,
  kDataPath2RdbRunIdType,
  kGamma1RdbRunIdType,
  kGamma2RdbRunIdType
};

/**
 * Parse a run_id string and return the formatted read-name prefix and the detected RDB run ID type.
 * @param run_id The run_id string from the RDB file header.
 * @return A tuple of (prefix string, RdbRunIdType).
 */
std::tuple<std::string, RdbRunIdType> BuildNamePrefix(const std::string& run_id);

/**
 * Encode a u32 cell index to a 6-character 0-padded base-64 string.
 *
 * Uses alphabet (0-9 A-Z a-z + _) with numeric base-64 conversion and
 * leading-zero padding so the output length is always exactly 6 characters.
 *
 * @param value The u32 cell index to encode.
 * @return 6-character base-64 string.
 */
std::string EncodeU32ToBase64(u32 value);

/**
 * The RdbSequenceReader is a \link SequenceReader and manages the process of
 * decoding RDB file content into batches. The RDB file is a binary file
 * containing sequencing output data. As this decoding involves moving back and
 * forth in the file (random access), it uses \link MemoryFileStream maintain
 * multiple memory cached file stream as needed for optimized reading of the
 * file.
 */
class RdbSequenceReader : public SequenceReader {
 public:
  explicit RdbSequenceReader(const fs::path& sequence_file_path);

  ~RdbSequenceReader() override = default;

  BatchStatistics ReadBatchIntoArena(FixedReadRecordBatch& batch) override;

 private:
  void LoadBlock();
  u32 ReadSingleFixed(FixedReadRecord& rec);
  bool HasMoreReads() const;

  std::ifstream _file;
  // the file path is stored here for logging purposes in case of errors during reading, so that we can log the file
  // path even if the error occurs in another thread
  fs::path _file_path;
  MemoryFileStream<> _read_seq_stream;
  MemoryFileStream<> _read_qual_stream;
  MemoryFileStream<> _read_len_stream;
  MemoryFileStream<> _read_cell_index_stream;

  u8 _qual_buf[kMaxReadLength]{};

  // current data block state
  std::streamsize _file_size_bytes;
  std::streamoff _cur_block_start;
  std::streamsize _cur_block_size;
  u32 _cur_batch_index;
  u32 _num_reads_in_block;
  u32 _cur_read_index;
  std::optional<u32> _cycle_id = std::nullopt;

  std::string _read_name_prefix;
};

}  // namespace xoos::demux
