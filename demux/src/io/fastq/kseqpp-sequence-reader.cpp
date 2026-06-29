#include "kseqpp-sequence-reader.h"

namespace xoos::demux {
KSeqppSequenceReader::KSeqppSequenceReader(const fs::path& path) : _seq_stream_in{path.c_str()} {}

BatchStatistics KSeqppSequenceReader::ReadBatchIntoArena(FixedReadRecordBatch& batch) {
  if (_seq_stream_in.eof()) {
    batch.num_records = 0;
    batch.num_bytes = 0;
    batch.num_bases = 0;
    return BatchStatistics{};
  }

  return _seq_stream_in.ReadBatchFixed(batch, batch.Capacity());
}
}  // namespace xoos::demux
