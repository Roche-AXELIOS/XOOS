#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <taskflow/taskflow.hpp>

#include "io/read-record.h"
#include "task/sink-worker.h"

namespace xoos::demux {

class FlowManager;

/**
 * Sink is responsible for writing the trimmed records to disk. It uses SinkWorker to write the records to disk.
 * The Sink class is thread-safe and can be used by multiple threads to write records to disk.
 *
 * Before demuxing starts, a sink is created for every possible SID; the constructor of the Datasink will create an
 * fastq output file for every SID - this means that a large number of files will be created that later turn out not
 * to contain any data; while this seems to be a waste of resources, it reduces the number of runtime memory allocations
 * which is always a good thing performance- and stability-wise.
 *
 * Multiple Formatting tasks can push FormattedRecords. To maximize throughput and minimize memory usage, there
 * are some ground rules:
 * 1) Assume that scalability is limited; at some point, throwing more CPUs at a problem will not increase throughput.
 *    With hundreds of input datasets, we therefore will limit the number of datasets that are processed in parallel.
 *    The number of simultaneous datasets processed is a function of thread::hardware_concurrency().
 * 2) This then means that we can limit the number of buffers that need to be allocated to accommodate the input data.
 *    We need at least one buffer per active SID per dataset (which usually amounts to roughly a dozen or so per
 *    dataset) that the formatting thread can write into, and a few buffers that the writer/compression threads can
 *    use to write the data out.
 * 3) We need to distinguish the two main cases: the writer thread(s) can keep up with the generated data, or they
 * can't. In the first case, we can keep the number of buffers to a minimum, and we can keep the memory usage low. In
 * the second case, adding more buffers to cache data will potentially cause unlimited memory usage, which is
 * unacceptable. In practice, we therefore either need to throttle the input data, or we need to add more writer threads
 * - but the total number of buffers that are used should be limited. With 1 MB buffers and typically only a few SID
 * active, we should strive to use a few dozen buffers at most - this will maximize throughput as data will stick around
 *    in L3 cache. Benchmarks should determine the sweet spot, but we'll start with 4 buffers per writer thread.
 * 4) We can allocate many more buffers upfront (say a 1000 buffers) to make sure we'll not run out of resources, but
 *    the key is that we need to limit the number of buffers that are actually used.
 *
 */

class Sink {
 public:
  Sink(const FlowManager& exec, SinkCompressionType compression_type, const std::filesystem::path& path,
       std::string_view input_name, std::string_view sample_name, size_t sid_nr,
       std::optional<size_t> compression_level, size_t max_number_of_workers);
  Sink(const Sink&) = delete;
  Sink& operator=(const Sink&) = delete;
  Sink(Sink&&) = delete;
  ~Sink() = default;

  void WriteData(const SinkData& data);

  size_t NumberOfWorkers() const { return _number_of_sink_workers; }

  size_t SidID() const { return _sid_nr; }

  size_t WriteTimeInNs() const { return _write_time; }

  tf::AsyncTask CreateWriteTask(const SinkData& data);
  void FlushAndClose();

 private:
  const FlowManager& _exec;
  std::vector<std::shared_ptr<SinkBase>> _sink_workers;
  const size_t _max_sink_workers;
  // When scheduling a new write task, we need to set the dependency of the previous write. Because only a few
  // writers are active, I decided to add the logic required for this to the Sink class.
  std::vector<tf::AsyncTask> _recent_tasks;
  std::atomic<size_t> _number_of_sink_workers{0};
  std::atomic<size_t> _push_count{0};
  std::atomic<size_t> _write_time{0};
  const std::filesystem::path& _path;
  const size_t _sid_nr;
  const u32 _compression_level;
  const SinkCompressionType _compression_type;
  const std::string _input_name;
  const std::string _sample_name;
  std::mutex _create_mutex;  // Prevents tasks from creating multiple output files at a time.

  friend class FlowManager;  // Allow FlowManager to call AddSinkWorker().
  void AddSinkWorker();
};

}  // namespace xoos::demux
