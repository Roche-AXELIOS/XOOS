#pragma once

#include <atomic>
#include <memory>
#include <vector>

#include "core/demux-and-trim-pipeline.h"
#include "io/sequence-reader.h"
#include "task/batch.h"

namespace xoos::demux {
class Sink;
class FlowManager;

/// @brief FlowContext holds the state of the flow graph associated with processing one input file.
class FlowContext {
 public:
  FlowContext(const FlowManager& mgr, const fs::path& input_file);
  ~FlowContext();

  /// Read data from the input file and write the data into the memory area pre-allocated for the batch.
  /// Return the number of sequences read (which usually equals the batch size)
  BatchStatistics GetBatch(size_t batch_nr);

  const FlowManager& GetManager() const { return _manager; }

  const DemuxAndTrimParam& DemuxParam() const { return _param; }

  size_t NumSequences() const { return _nr_sequences; }

  size_t NumBases() const { return _nr_bases; }

  size_t NumBytes() const { return _nr_bytes; }

  std::chrono::nanoseconds::rep ReadTimeInNs() const { return _read_time; }

  void ReportStatistics() const;

  FixedReadRecordBatch& GetBatchData(size_t batch_nr) { return *_batch_data[BatchNrToIndex(batch_nr)]; }

  size_t AddToDemuxTime(uint64_t time) { return _demux_time += time; }

  size_t AddToFormatTime(uint64_t time) { return _format_time += time; }

  size_t AddToAlignmentTime(uint64_t time) { return _alignment_time += time; }

  // Used by the write task to signal that a batch has been completed and was handed over to a sink.
  void MarkBatchCompleted() { ScheduleNewBatch(); }

  // Set writing task for a batch number.
  void SetWriteTask(size_t batch_nr, tf::AsyncTask write_task);

  // Write output data associated with a batch number.
  void WriteData(const SinkData& sink_data);

  // Get the write task that was used to write an earlier batch. This is used to chain write tasks to write
  // out data such that only one thread accesses the writer at a time.
  tf::AsyncTask CreateWriteTask(const SinkData& data) const;

  // For determining when to stop processing
  std::atomic<size_t> nr_writes_scheduled{0ul};
  std::atomic<size_t> nr_tasks_scheduled{0ul};

 private:
  const FlowManager& _manager;
  const DemuxAndTrimParam& _param;
  std::vector<std::shared_ptr<BatchTask>> _tasks;
  /// @brief Memory used to process the batches; every active batch needs its own memory area.
  std::vector<std::unique_ptr<FixedReadRecordBatch>> _batch_data;
  std::vector<std::unique_ptr<Sink>> _sinks;  // One sink per valid SID
  std::shared_ptr<SequenceReader> _reader;
  const fs::path& _input_file;
  std::atomic<size_t> _next_batch_id{0};
  std::atomic<std::chrono::nanoseconds::rep> _read_time{0ul};
  std::atomic<uint64_t> _demux_time{0ul};
  std::atomic<uint64_t> _alignment_time{0ul};
  std::atomic<uint64_t> _format_time{0ul};
  std::atomic<bool> _input_available{true};
  std::atomic<bool> _first_batches_scheduled{false};

  // For tracking the number of sequences read from the input file.
  std::atomic<size_t> _nr_sequences{0ul};
  std::atomic<size_t> _nr_bases{0ul};
  const size_t _nr_bytes;

  // Required to determine when to wrap up processing
  std::atomic<size_t> _nr_writes_completed{0ul};

  /// This is one of the main functions of the FlowContext; it schedules a future batch to be processed. More
  /// details are in the implementation file.
  void ScheduleNewBatch();

  /// Schedule the first batches to be processed; these need special treatment because not all dependencies are
  /// available yet. More details are in the implementation file.
  void ScheduleFirstBatches();

  // Helper function to convert a batch number to an index in the tasks vector.
  size_t static BatchNrToIndex(size_t batch_nr) { return batch_nr % BatchTask::kCircularBufferSize; }
};
}  // namespace xoos::demux
