#include "task/sink.h"

#include <fmt/format.h>
#include <xoos/error/error.h>

#include "task/batch.h"
#include "task/flow-manager.h"

namespace xoos::demux {
Sink::Sink(const FlowManager& exec, const SinkCompressionType compression_type, const std::filesystem::path& path,
           const std::string_view input_name, const std::string_view sample_name, const size_t sid_nr,
           const std::optional<size_t> compression_level, const size_t max_number_of_workers)
    : _exec(exec),
      _sink_workers(BatchTask::kMaxNumberOfSinkWorkers, nullptr),
      _max_sink_workers(max_number_of_workers),
      _recent_tasks(BatchTask::kMaxNumberOfSinkWorkers),
      _path(path),
      _sid_nr(sid_nr),
      _compression_level([&compression_type, &compression_level]() -> u32 {
        if (compression_level.has_value()) {
          const auto level = compression_level.value();
          if (compression_type == SinkCompressionType::kGzip) {
            if (level < kMinCompressionLevel || level > kMaxGzipCompressionLevel) {
              throw error::Error("Invalid gzip compression level: {}, must be in the range [{}, {}]", level,
                                 kMinCompressionLevel, kMaxGzipCompressionLevel);
            }
          } else if (compression_type == SinkCompressionType::kZstd &&
                     (level < kMinCompressionLevel || level > kMaxZstdCompressionLevel)) {
            throw error::Error("Invalid zstd compression level: {}, must be in the range [{}, {}]", level,
                               kMinCompressionLevel, kMaxZstdCompressionLevel);
          }
          return static_cast<u32>(level);
        }
        // No value provided, use defaults
        if (compression_type == SinkCompressionType::kGzip) {
          return kDefaultGZipCompressionLevel;
        }
        if (compression_type == SinkCompressionType::kZstd) {
          return kDefaultZstdCompressionLevel;
        }
        // For kNone or other types
        return 0;
      }()),
      _compression_type(compression_type),
      _input_name(input_name),
      _sample_name(sample_name) {}

void Sink::WriteData(const SinkData& data) {
  auto sink_worker = _sink_workers[data.batch_id % _number_of_sink_workers];
  if (data.length) {
    _write_time += sink_worker->WriteData(data);
    ++_push_count;
  }
}

void Sink::FlushAndClose() {
  for (auto& sink_worker : _sink_workers) {
    if (sink_worker) {
      _write_time += sink_worker->FlushAndClose();
    }
  }
}

tf::AsyncTask Sink::CreateWriteTask(const SinkData& data) {
  // fancy scaling stuff happens here later
  if (_number_of_sink_workers == 0) {
    AddSinkWorker();

  } else {
    if (_push_count > 1000 && _number_of_sink_workers == 1) {
      for (size_t i = 1ul; i < BatchTask::kMaxNumberOfSinkWorkers; ++i) {
        AddSinkWorker();
      }
    }
  }

  auto& recent_task = _recent_tasks[data.batch_id % _number_of_sink_workers];
  if (recent_task.empty() || recent_task.is_done()) {
    // Previous task is done or not present, no dependency needed
    recent_task = _exec.Executor().silent_dependent_async(data);
  } else {
    recent_task = _exec.Executor().silent_dependent_async(data, recent_task);
  }
  return recent_task;
}

void Sink::AddSinkWorker() {
  std::scoped_lock lock(_create_mutex);
  // Only create a new sink worker if we have not yet reached the configured maximum.
  // This bounds the number of writer threads/output files per sink to _max_sink_workers.
  if (_number_of_sink_workers < _max_sink_workers) {
    const auto sink_nr{static_cast<size_t>(_number_of_sink_workers)};

    const auto extension = _compression_type == SinkCompressionType::kZstd ? ".fastq.zst"
                           : (_compression_level == 0 && _compression_type == SinkCompressionType::kGzip) ||
                                   _compression_type == SinkCompressionType::kNone
                               ? ".fastq"
                               : ".fastq.gz";

    auto output_path = _path / _sample_name;
    // we only add sink nr if we have multiple workers
    if (_max_sink_workers > 1) {
      output_path = output_path / fmt::format("{}.{:02}{}", _input_name, sink_nr, extension);
    } else {
      output_path = output_path / fmt::format("{}{}", _input_name, extension);
    }

    if (_sid_nr == 0) {
      // special case for unassigned reads
      output_path = _path / fmt::format("{}.{:02}{}", _input_name, sink_nr, extension);
    }
    if (fs::exists(output_path)) {
      if (_exec.param.overwrite) {
        // delete the file if it already exists and we are allowed to overwrite
        fs::remove(output_path);
      } else {
        throw error::Error(
            fmt::format("Output file {} already exists. Specify --overwrite to allow overwriting existing files.",
                        output_path.string()));
      }
    }

    // create a directory for the output_path parent if it does not exist
    std::filesystem::create_directories(output_path.parent_path());

    switch (_compression_type) {
      case SinkCompressionType::kNone:
        _sink_workers[_number_of_sink_workers] = std::make_shared<SinkRaw>(output_path);
        break;
      case SinkCompressionType::kGzip:
        _sink_workers[_number_of_sink_workers] = std::make_shared<SinkGZip>(output_path, _compression_level);
        break;
      case SinkCompressionType::kZstd:
        _sink_workers[_number_of_sink_workers] = std::make_shared<SinkZStd>(output_path, _compression_level);
        break;
      default:  // Null device
        _sink_workers[_number_of_sink_workers] = std::make_shared<SinkNull>(output_path);
        break;
    }
    ++_number_of_sink_workers;
    // Copy over most recently used task.
    if (_number_of_sink_workers > 1) {
      _recent_tasks[_number_of_sink_workers - 1] = _recent_tasks[_number_of_sink_workers - 2];
    }
  }
}
}  // namespace xoos::demux
