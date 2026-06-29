#pragma once

#include <taskflow/taskflow.hpp>

#include "io/read-record.h"
#include "task/alignment.h"
#include "task/demux.h"
#include "task/formatter.h"
#include "task/reader.h"

namespace xoos::demux {

/// BatchTask is a struct that contains all the tasks needed to process a batch of data; in addition to the state
/// information of reader, demux, and formatter, it also contains the four taskflow tasks that take
/// care of running the tasks when the time comes.
/// We cannot schedule the writing of output before we have finished processing the first
/// three tasks, because we need to know which SIDs were encountered. After Formatting the data, we can schedule
/// writing output data to the correct Sink (one output file per SID per input dataset) and we can do this in parallel.
/// These sinks are created by the Formatter object, and all sinks are run in parallel before we can wrap up the batch.
/// The writer_task is a dummy task that depends on the completion of the sink; it is not actually doing any work, but
/// used for synchronization so that we can be sure that all the sinks have been created before we finish the batch.

struct BatchTask {
  /// Important constants that control the scheduler behavior and the memory usage:
  /// - kBatchesAhead: The number of batches to be scheduled ahead of the current batch (the read task of batch N
  ///   schedules batch N + kBatchesAhead).
  /// - kBatchSinkDelay: The time (expressed in batches) that we allow the sink threads to write out the data. Of
  ///   course, other threads will be doing work in the meantime.
  /// - kMaxNumberOfSinkWorkers: Maximum parallelism for writing out the data. While this number does not
  ///   directly affect scheduling, it does affect data dependencies in the circular buffer.
  /// - kCircularBufferSize: The size of the circular buffer that holds all data.

  constexpr static size_t kBatchesAhead = 10;
  constexpr static size_t kBatchSinkDelay = 4;
  constexpr static size_t kMaxNumberOfSinkWorkers = 8;
  constexpr static size_t kCircularBufferSize = kBatchesAhead      // Task N schedules task N + kBatchesAhead
                                                + kBatchesAhead    // Task N depends on format of task N - kBatchesAhead
                                                + kBatchSinkDelay  // account for sink latency
                                                + kMaxNumberOfSinkWorkers;  // account for sink parallelism

  /// Constructor that initializes the reader, demux, and formatter objects with the given batch number.
  BatchTask(FlowContext& exec, size_t batch_nr)
      : reader(exec, batch_nr), demux(exec, batch_nr), alignment(exec, batch_nr), formatter(exec, batch_nr) {}

  Reader reader;                // State information and operator() for reading data.
  Demux demux;                  // State information and operator() for demultiplexing data.
  PairwiseAlignment alignment;  // State information and operator() for aligning data (only for HD Demux)
  Formatter formatter;  // State information and operator() for formatting data; this includes creating the sinks!

  tf::AsyncTask reader_task;      // Taskflow task for reading data.
  tf::AsyncTask demux_task;       // Taskflow task for demultiplexing data.
  tf::AsyncTask alignment_task;   // Taskflow task for aligning data (only for HD Demux)
  tf::AsyncTask formatting_task;  // Taskflow task for formatting data.
  tf::AsyncTask writer_task;      // Dummy taskflow task; depends on completion of all sink tasks created by formatter.
};

}  // namespace xoos::demux
