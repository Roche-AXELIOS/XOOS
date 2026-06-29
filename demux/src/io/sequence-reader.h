#pragma once

#include <filesystem>
#include <memory>

#include "read-record.h"

namespace fs = std::filesystem;

namespace xoos::demux {
struct BatchStatistics {
  size_t num_sequences{0UL};
  size_t num_bases{0UL};
};

class SequenceReader {
 public:
  SequenceReader() = default;
  virtual ~SequenceReader() = default;

  /// Reads the next batch of records into the provided arena. It is assumed that the batch memory
  /// already has been pre-allocated and the desired number of records (i.e. batch size) is determined
  /// by the number of pre-allocated records in batch. This function updates the contents of batch and
  /// returns the number of records read (which is usually equal to the number of pre-allocated records).
  virtual BatchStatistics ReadBatchIntoArena(FixedReadRecordBatch& batch) = 0;
};

bool IsSequenceFileFormat(const fs::path& sequence_file_path);

std::shared_ptr<SequenceReader> OpenSequenceFile(const fs::path& sequence_file_path);

}  // namespace xoos::demux
