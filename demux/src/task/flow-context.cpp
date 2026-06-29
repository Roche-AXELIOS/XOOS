#include "task/flow-context.h"

#include <fmt/format.h>
#include <xoos/log/logging.h>

#include <algorithm>
#include <utility>

#include "task/flow-manager.h"

namespace xoos::demux {
FlowContext::FlowContext(const FlowManager& mgr, const fs::path& input_file)
    // max_active_batches is the number of batches that are processed in parallel; with more batches, the scheduler
    // can optimize the task graph better at the cost of memory consumption.
    : _manager(mgr),
      _param(mgr.param),
      _input_file(input_file),
      _next_batch_id(BatchTask::kBatchesAhead),
      _nr_bytes(file_size(input_file)) {
  _tasks.reserve(BatchTask::kCircularBufferSize);

  // Allocate memory for the tasks; it needs to be allocated before we start scheduling the tasks.
  for (size_t batch_nr = 0; batch_nr < BatchTask::kCircularBufferSize; ++batch_nr) {
    _batch_data.emplace_back(std::make_unique<FixedReadRecordBatch>(_manager.param.batch_size));
  }

  // Allocate a potential write buffer for every valid SID. Assume the minimum valid SID ID is 1, add an extra
  // buffer at position 0 for unmapped reads.
  uint max_id{0};
  for (const auto& [sid_id, sid] : mgr.sid_pool) {
    max_id = std::max(max_id, sid_id);
  }

  // 0 is reserved for unknown SIDs
  _sinks.resize(max_id + 2);
  for (const auto& [sid_id, sid] : mgr.sid_pool) {
    // Create an output sink for every valid SID. Note that the SID ID is 1-based.
    // To create a unique output filename, take the stem of the input file (so we have a hunch which input file was
    // used + a hash value of the input file to make sure that the output files are unique.
    // To limit string lengths, we limit the hash value length.
    std::string prefix = fmt::format("{}.{}_{}", sid.name, input_file.stem().string(),
                                     (std::hash<std::string>{}(input_file.string())) % 100000000UL);

    _sinks[sid_id + 1] =
        std::make_unique<Sink>(_manager, _param.compression_type, _param.out_dir, prefix, sid.name, sid_id + 1,
                               _param.compression_level, _param.writing_threads_per_sample);
  }

  // Dummy sink for all failed reads (for any reason).
  std::string prefix = fmt::format("raw_failed.{}_{}", input_file.stem().string(),
                                   std::hash<std::string>{}(input_file.string()) % 100000000UL);
  _sinks[0] = std::make_unique<Sink>(_manager, _param.compression_type, _param.out_dir, prefix, "raw_failed", 0,
                                     _param.compression_level, _param.writing_threads_per_sample);

  try {
    _reader = OpenSequenceFile(input_file);
  } catch (const std::exception& e) {
    Logging::Error("Error opening input file {}: {}", input_file.string(), e.what());
    Task::SetTaskException(std::current_exception());
    return;
  }

  ScheduleFirstBatches();
}

tf::AsyncTask FlowContext::CreateWriteTask(const SinkData& data) const {
  // Forward this request to the sink that will be used.
  return _sinks[data.sid_id]->CreateWriteTask(data);
}

void SinkData::operator()() const { flow_context->WriteData(*this); }

void FlowContext::WriteData(const SinkData& data) {
  auto& sink = _sinks[data.sid_id];
  sink->WriteData(data);
  ++_nr_writes_completed;
}

FlowContext::~FlowContext() {  // don't expect this to happen often
  while (!Task::HasException() && (_input_available || _nr_writes_completed != nr_writes_scheduled)) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  if (Task::HasException()) {
    // If there is an exception we want to exit without writing out anything
    // and not throwing an error in the constructor
    return;
  }

  for (auto& s : _sinks) {
    if (s) {
      s->FlushAndClose();
    }
  }
  ReportStatistics();
}

// This is perhaps the most important function in the whole program; it schedules new tasks for a future batch. When
// figuring out all the cross-dependencies between tasks, I created a high-level overview of the task graph using
// GraphViz. This helped me to understand the dependencies between tasks and to make sure that the tasks were scheduled
// in the right order. Here is an example of one of the tasks in the GraphViz input:
//  subgraph cluster_16 {
//    "read N+16"->"demux N+16"->"format N+16"
//    "read N+15"->"read N+16"
//    "write N+8"->"read N+16"
//    "format N+12"->"read N+16"
//    label = "batch N+16";
//  }
// The subgraph contains the tasks for batch N+16. The arrows indicate the dependencies between the tasks. The read task
// is dependent on the read task of the previous batch (read N+15), the write task of the oldest batch (write N+8), and
// the format task of the previous batch (format N+12).
// This function creates all the dependencies between the tasks for a future batch; the format task of the N+12 task
// will generate the sink tasks and writer object that manage the output files.
//
void FlowContext::ScheduleNewBatch() {
  if (!_input_available) {
    return;  // stop scheduling new batches if we have reached the end of the input file.
  }

  while (!_first_batches_scheduled.load()) {
    // empty body
  }

  // Schedule task for the next "future" batch. The new batch number is the current value of next_batch_id.
  // Please realize that the reader that calls this function has itself the batch number batch_nr - nr_batches_ahead.

  // Let's consider the edge case where we have already scheduled nr_batches_ahead batches. That means that the
  // next batch ID that would be scheduled is nr_batches_ahead (e.g., 4); the batch we're about to schedule should
  // have a dependency on the Formatter task of batch_nr - nr_batches_ahead (e.g. 0).
  size_t batch_nr = _next_batch_id++;  // does std::atomic_fetch_add under the hood

  // add dependency on previous reader and the formatter of the previous "oldest" batch
  auto previous_read_index{BatchNrToIndex(batch_nr - 1)};
  auto previous_format_index{BatchNrToIndex(batch_nr - BatchTask::kBatchesAhead)};

  // Create a new batch task and add dependencies on the previous batches.
  auto p_task{std::make_shared<BatchTask>(*this, batch_nr)};

  // Overwrite the oldest task with the new one - the BatchNrToIndex function will wrap around and basically create
  // a circular buffer. This is a simple way to keep the number of active tasks in memory constant.
  _tasks[BatchNrToIndex(batch_nr)] = p_task;

  auto& executor = _manager.Executor();

  // The reader task is dependent on the reader task of the previous batch, the format task of batch #
  // batch_nr - BatchTask::kBatchesAhead (which is guaranteed to be available).
  // I also use dependencies on the writer task of a previous batch; I do
  // realize that that dependency is already implicit because I also set up dependencies in the sink tasks that
  // are created, but the additional dependency does not seem to hurt performance - actually, performance seems
  // to be a tad better because the pipeline becomes "more predictable".

  constexpr auto kWriterDelay = BatchTask::kBatchesAhead + BatchTask::kBatchSinkDelay;
  auto previous_writer_batch_nr = batch_nr < kWriterDelay ? 0ul : batch_nr - kWriterDelay;
  auto previous_writer_index{BatchNrToIndex(previous_writer_batch_nr)};

  // We can now initialize the various tasks using the correct dependencies.
  p_task->reader_task = executor.silent_dependent_async(
      p_task->reader.Name(), p_task->reader, _tasks[previous_read_index]->reader_task,
      _tasks[previous_format_index]->formatting_task, _tasks[previous_writer_index]->writer_task);

  // Demuxing is dependent on the reader task.
  p_task->demux_task = executor.silent_dependent_async(p_task->demux.Name(), p_task->demux, p_task->reader_task);

  // Alignment is dependent on the demuxing task - but only for HD Duplex data.
  if (_manager.is_duplex) {
    p_task->alignment_task =
        executor.silent_dependent_async(p_task->alignment.Name(), p_task->alignment, p_task->demux_task);
    // align task is scheduled.
    nr_tasks_scheduled += 1;
  }

  // Formatting is dependent on the demuxing task. We allow for out-of-order writes, so we do not have to create a
  // dependency on _tasks[previous_read_index]->formatting_task.
  p_task->formatting_task = executor.silent_dependent_async(
      p_task->formatter.Name(), p_task->formatter, _manager.is_duplex ? p_task->alignment_task : p_task->demux_task,
      _tasks[previous_read_index]->formatting_task);

  // reader, demux, and formatter tasks are scheduled.
  nr_tasks_scheduled += 3;

  // All tasks for the new batch have now been scheduled. The dependencies associated with the writer task are
  // generated during the formatting task, we do not know beforehand which SIDs are present in each batch.
  // See formatting task for more details about how the sink tasks and writer task are created.
}

void FlowContext::ScheduleFirstBatches() {
  // We need to schedule the first few batches manually.

  // Start our blocking task for the first few batches.
  auto& executor = _manager.Executor();

  std::shared_ptr<BatchTask> prev_task;
  for (size_t batch_nr = 0; batch_nr < BatchTask::kBatchesAhead; ++batch_nr) {
    auto p_task{std::make_shared<BatchTask>(*this, batch_nr)};
    _tasks.emplace_back(p_task);
    // Block the read tasks from executing until we unblock; this also means that any dependent tasks will not execute.
    // So: we only need to block the first reader task as subsequent ones depend on it.
    if (batch_nr == 0) {  // no dependency at all except for blocker.
      p_task->reader_task = executor.silent_dependent_async(p_task->reader.Name(), p_task->reader);
    } else {
      p_task->reader_task =
          executor.silent_dependent_async(p_task->reader.Name(), p_task->reader, prev_task->reader_task);
    }

    p_task->demux_task = executor.silent_dependent_async(p_task->demux.Name(), p_task->demux, p_task->reader_task);

    if (_manager.is_duplex) {
      p_task->alignment_task =
          executor.silent_dependent_async(p_task->alignment.Name(), p_task->alignment, p_task->demux_task);
      // align task is scheduled.
      nr_tasks_scheduled += 1;
    }

    auto previous_task{_manager.is_duplex ? p_task->alignment_task : p_task->demux_task};

    if (batch_nr == 0) {
      p_task->formatting_task =
          executor.silent_dependent_async(p_task->formatter.Name(), p_task->formatter, previous_task);
    } else {
      p_task->formatting_task = executor.silent_dependent_async(p_task->formatter.Name(), p_task->formatter,
                                                                previous_task, prev_task->formatting_task);
    }

    // reader, demux, and formatter tasks are scheduled.
    nr_tasks_scheduled += 3;
    prev_task = p_task;
  }

  // We also will insert dummy tasks for all other tasks in the circular buffer; these will be overwritten by
  // subsequent batches.
  for (size_t batch_nr = BatchTask::kBatchesAhead; batch_nr < BatchTask::kCircularBufferSize; ++batch_nr) {
    _tasks.emplace_back(std::make_shared<BatchTask>(*this, batch_nr));
  }
  _first_batches_scheduled.store(true);
}

void FlowContext::SetWriteTask(const size_t batch_nr, tf::AsyncTask write_task) {
  _tasks[BatchNrToIndex(batch_nr)]->writer_task = std::move(write_task);
}

// Read data from the input file.
BatchStatistics FlowContext::GetBatch(const size_t batch_nr) {
  try {
    // Called by the reader task to read a batch of data from the input file.
    auto& buffer{*_batch_data[BatchNrToIndex(batch_nr)]};
    auto batch_stats{_reader->ReadBatchIntoArena(buffer)};
    _input_available = batch_stats.num_sequences > 0;
    _nr_bases += batch_stats.num_bases;
    _nr_sequences += batch_stats.num_sequences;
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(buffer.end_time - buffer.start_time);
    _read_time += duration.count();
    return batch_stats;
  } catch (const std::exception& e) {
    _input_available = false;
    throw;
  }
}

void FlowContext::ReportStatistics() const {
  size_t write_time{0};
  for (auto& sink : _sinks) {
    if (sink) {
      write_time += sink->WriteTimeInNs();
    }
  }

  auto total_tasks = nr_tasks_scheduled + nr_writes_scheduled;

  if (_manager.is_duplex) {
    Logging::Info(
        "{} ({} bytes): {} reads, {} bases ==> read: {:.2f} ms, demux: {:.2f} ms, alignment: {:.2f} ms, format: {:.2f} "
        "ms, write: {:.2f} ms "
        "({} tasks)",
        _input_file.string(), _nr_bytes, static_cast<size_t>(_nr_sequences), static_cast<size_t>(_nr_bases),
        (1e-6f * static_cast<f32>(_read_time)), (1e-6f * static_cast<f32>(_demux_time)),
        (1e-6f * static_cast<f32>(_alignment_time)), (1e-6f * static_cast<f32>(_format_time)),
        (1e-6f * static_cast<f32>(write_time)), total_tasks);
  } else {
    Logging::Info(
        "{} ({} bytes): {} reads, {} bases ==> read: {:.2f} ms, demux: {:.2f} ms, format: {:.2f} ms, write: {:.2f} ms "
        "({} tasks)",
        _input_file.string(), _nr_bytes, static_cast<size_t>(_nr_sequences), static_cast<size_t>(_nr_bases),
        (1e-6f * static_cast<f32>(_read_time)), (1e-6f * static_cast<f32>(_demux_time)),
        (1e-6f * static_cast<f32>(_format_time)), (1e-6f * static_cast<f32>(write_time)), total_tasks);
  }
}

}  // namespace xoos::demux
