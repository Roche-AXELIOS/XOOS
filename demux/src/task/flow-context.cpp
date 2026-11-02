#include "task/flow-context.h"

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include "task/flow-manager.h"

namespace xoos::demux {
FlowContext::FlowContext(const FlowManager& mgr, const fs::path& input_file, size_t input_index)
    // max_active_batches is the number of batches that are processed in parallel; with more batches, the scheduler
    // can optimize the task graph better at the cost of memory consumption.
    : _manager(mgr),
      _param(mgr.Param()),
      _input_file(input_file),
      _next_batch_id(BatchTask::kBatchesAhead),
      _nr_bytes(file_size(input_file)) {
  _tasks.reserve(BatchTask::kCircularBufferSize);

  // Allocate memory for the tasks; it needs to be allocated before we start scheduling the tasks.
  // NOSONAR(M23_002) -- loop over a positive compile-time constant bound; the first comparison being true is
  // expected loop behaviour, not a redundant condition.
  for (size_t batch_nr = 0; batch_nr < BatchTask::kCircularBufferSize; ++batch_nr) {  // NOSONAR(M23_002)
    const auto budget = HasTrimmingDebugMode(_manager.Param().trimming_debug, TrimmingDebugMode::kRawSeq)
                            ? kBufferSize + kMaxReadLength
                            : kBufferSize;
    _batch_data.emplace_back(std::make_unique<FixedReadRecordBatch>(_manager.Param().batch_size, budget));
  }

  // Allocate a sink buffer for every SID plus one extra at index 0 for unassigned/failed reads.
  // SID IDs are 0-based; sink index is sid.id + 1.
  // In kAllSplit mode, an additional set of sinks is created for partial reads at offset sid_pool.size().
  if (mgr.SidPoolValue().empty()) [[unlikely]] {
    throw error::Error("sid_pool must contain at least one SID");
  }
  const bool split_mode = _param.read_length_mode == ReadLengthMode::kAllSplit;
  const auto sinks_per_category = mgr.SidPoolValue().size();
  _sinks.resize(1 + sinks_per_category + (split_mode ? sinks_per_category : 0));
  // use full path for input name to prevent collisions
  const auto input_filename = input_file.string();

  for (const auto& sid : mgr.SidPoolValue()) {
    // The prefix includes the input file's index in the sorted input list to guarantee unique temporary filenames
    // across input files (even when multiple files share the same stem, e.g., files in different directories).
    // These temporary names are replaced with Everest-format names (<sample_name>-<part_index>) during the
    // post-processing rename step.
    std::string prefix = fmt::format("{}.{}_{}", sid.name, input_file.stem().string(), input_index);
    // Track the mapping from output file prefix to input file path for post-processing rename and provenance.
    mgr.RegisterOutputPrefix(prefix, input_filename);

    // In split mode, the normal sink index holds full reads under <sample>/full/.
    // In non-split modes, the sink writes directly to <sample>/.
    const auto sample_dir = split_mode ? fmt::format("{}/{}", sid.name, kFullReadSubdir) : sid.name;
    _sinks[sid.id + 1] = std::make_unique<Sink>(_manager,
                                                _param.compression_type,
                                                _param.out_dir,
                                                prefix,
                                                sample_dir,
                                                sid.id + 1,
                                                _param.compression_level,
                                                _param.writing_threads_per_sample);

    if (split_mode) {
      // Partial-read sinks at offset sinks_per_category, writing to <sample>/partial/.
      const auto partial_sample_dir = fmt::format("{}/{}", sid.name, kPartialReadSubdir);
      const auto partial_sink_index = sid.id + 1 + sinks_per_category;
      _sinks[partial_sink_index] = std::make_unique<Sink>(_manager,
                                                          _param.compression_type,
                                                          _param.out_dir,
                                                          prefix,
                                                          partial_sample_dir,
                                                          partial_sink_index,
                                                          _param.compression_level,
                                                          _param.writing_threads_per_sample);
    }
  }

  // Dummy sink for all failed reads (for any reason).
  std::string prefix = fmt::format("{}.{}_{}", kReservedRawFailedName, input_file.stem().string(), input_index);
  // Track the mapping from output file prefix to input file path for post-processing rename and provenance.
  mgr.RegisterOutputPrefix(prefix, input_filename);
  _sinks[0] = std::make_unique<Sink>(_manager,
                                     _param.compression_type,
                                     _param.out_dir,
                                     prefix,
                                     kReservedRawFailedName,
                                     0,
                                     _param.compression_level,
                                     _param.writing_threads_per_sample);

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

void SinkData::operator()() const {
  try {
    flow_context->WriteData(*this);
  } catch (const std::exception& e) {  // intentionally generic to prevent stall on any write failure
    // only log the first exception to avoid spamming the logs
    if (!Task::HasException()) [[unlikely]] {  // NOLINT(readability/braces)
      Logging::Error("Write task failed for SID {}: {}", sid_id, e.what());
      Task::SetTaskException(std::current_exception());
    }
    // Increment so the FlowContext destructor spin-wait can terminate.
    // Without this, _nr_writes_completed never reaches nr_writes_scheduled
    // and the destructor loops forever.
    flow_context->IncrementWritesCompleted();
  } catch (...) {  // catch all to prevent stall on any write failure, even non-standard exceptions
    if (!Task::HasException()) [[unlikely]] {  // NOLINT(readability/braces)
      Logging::Error("Write task failed for SID {}: unknown exception", sid_id);
      Task::SetTaskException(std::current_exception());
    }
    flow_context->IncrementWritesCompleted();
  }
}

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

  // Empty-sample materialization is handled globally in FlowManager::MaterializeEmptySamples(), not here: a FlowContext
  // is per input file, so doing it here would emit one empty file per input for a globally-empty sample.
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
  const auto previous_writer_batch_nr = batch_nr < kWriterDelay ? 0UL : batch_nr - kWriterDelay;
  auto previous_writer_index{BatchNrToIndex(previous_writer_batch_nr)};

  // We can now initialize the various tasks using the correct dependencies.
  p_task->reader_task = executor.silent_dependent_async(p_task->reader.Name(),
                                                        p_task->reader,
                                                        _tasks[previous_read_index]->reader_task,
                                                        _tasks[previous_format_index]->formatting_task,
                                                        _tasks[previous_writer_index]->writer_task);

  // Demuxing is dependent on the reader task.
  p_task->demux_task = executor.silent_dependent_async(p_task->demux.Name(), p_task->demux, p_task->reader_task);

  // Alignment is dependent on the demuxing task - but only for HD Duplex data.
  if (_manager.IsDuplex()) {
    p_task->alignment_task =
        executor.silent_dependent_async(p_task->alignment.Name(), p_task->alignment, p_task->demux_task);
    // align task is scheduled.
    nr_tasks_scheduled += 1;
  }

  // Formatting is dependent on the demuxing task. We allow for out-of-order writes, so we do not have to create a
  // dependency on _tasks[previous_read_index]->formatting_task.
  p_task->formatting_task =
      executor.silent_dependent_async(p_task->formatter.Name(),
                                      p_task->formatter,
                                      _manager.IsDuplex() ? p_task->alignment_task : p_task->demux_task,
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
  // NOSONAR(M23_002) -- loop over a positive compile-time constant bound; the first comparison being true is
  // expected loop behaviour, not a redundant condition.
  for (size_t batch_nr = 0; batch_nr < BatchTask::kBatchesAhead; ++batch_nr) {  // NOSONAR(M23_002)
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

    if (_manager.IsDuplex()) {
      p_task->alignment_task =
          executor.silent_dependent_async(p_task->alignment.Name(), p_task->alignment, p_task->demux_task);
      // align task is scheduled.
      nr_tasks_scheduled += 1;
    }

    auto previous_task{_manager.IsDuplex() ? p_task->alignment_task : p_task->demux_task};

    if (batch_nr == 0) {
      p_task->formatting_task =
          executor.silent_dependent_async(p_task->formatter.Name(), p_task->formatter, previous_task);
    } else {
      p_task->formatting_task = executor.silent_dependent_async(
          p_task->formatter.Name(), p_task->formatter, previous_task, prev_task->formatting_task);
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
  for (const auto& sink : _sinks) {
    if (sink) {
      write_time += sink->WriteTimeInNs();
    }
  }

  auto total_tasks = nr_tasks_scheduled + nr_writes_scheduled;

  if (_manager.IsDuplex()) {
    Logging::InfoNoExcept(
        "{} ({} bytes): {} reads, {} bases ==> read: {:.2f} ms, demux: {:.2f} ms, alignment: {:.2f} ms, format: {:.2f} "
        "ms, write: {:.2f} ms "
        "({} tasks)",
        _input_file.native(),
        _nr_bytes,
        static_cast<size_t>(_nr_sequences),
        static_cast<size_t>(_nr_bases),
        (1e-6F * static_cast<f32>(_read_time)),
        (1e-6F * static_cast<f32>(_demux_time)),
        (1e-6F * static_cast<f32>(_alignment_time)),
        (1e-6F * static_cast<f32>(_format_time)),
        (1e-6F * static_cast<f32>(write_time)),
        total_tasks);
  } else {
    Logging::InfoNoExcept(
        "{} ({} bytes): {} reads, {} bases ==> read: {:.2f} ms, demux: {:.2f} ms, format: {:.2f} ms, write: {:.2f} ms "
        "({} tasks)",
        _input_file.native(),
        _nr_bytes,
        static_cast<size_t>(_nr_sequences),
        static_cast<size_t>(_nr_bases),
        (1e-6F * static_cast<f32>(_read_time)),
        (1e-6F * static_cast<f32>(_demux_time)),
        (1e-6F * static_cast<f32>(_format_time)),
        (1e-6F * static_cast<f32>(write_time)),
        total_tasks);
  }
}

}  // namespace xoos::demux
