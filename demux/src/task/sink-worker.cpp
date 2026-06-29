#include "task/sink-worker.h"

#include <fmt/format.h>
#include <xoos/error/error.h>

#include <fstream>

#include "utility/stop-watch.h"

namespace xoos::demux {

SinkRaw::SinkRaw(const std::filesystem::path& path) : SinkBase(path) {
  _output_file_stream.open(path, std::ios::binary);
  if (!_output_file_stream.is_open()) {
    throw error::Error("Failed to open file: {}", path.string());
  }
}

SinkRaw::~SinkRaw() { FlushAndClose(); }

size_t SinkRaw::WriteData(const SinkData& data) {
  StopWatch sw;
  _output_file_stream.write(data.p_data, static_cast<std::streamsize>(data.length));
  if (!_output_file_stream.good()) {
    throw error::Error("Failed to write to file: {}", output_file_name);
  }
  _output_file_stream.flush();
  return sw.ElapsedTime();
}

size_t SinkRaw::FlushAndClose() {
  auto start{std::chrono::high_resolution_clock::now()};
  if (_output_file_stream.is_open()) {
    _output_file_stream.flush();
    _output_file_stream.close();
  }
  auto end{std::chrono::high_resolution_clock::now()};
  return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

SinkGZip::SinkGZip(const std::filesystem::path& path, size_t compression_level) : SinkBase(path) {
  const std::string write_mode = fmt::format("w{}", compression_level);
  _output_file = gzopen(path.c_str(), write_mode.c_str());
  if (_output_file == Z_NULL) {
    throw error::Error("Failed to open file: '{}'", path.string());
  }
}

SinkGZip::~SinkGZip() { FlushAndClose(); }

size_t SinkGZip::WriteData(const SinkData& data) {
  auto start{std::chrono::high_resolution_clock::now()};

  size_t nr_written = gzwrite(_output_file, data.p_data, data.length);
  if (nr_written != data.length) {
    int err_num;
    const char* error_string = gzerror(_output_file, &err_num);
    if (err_num == Z_ERRNO) {
      throw error::Error("Failed to write to file: {}. System error: {}", output_file_name, error_string);
    }
    throw error::Error("gzip error: {}. Failed to write to file: {}", std::string(error_string), output_file_name);
  }
  gzflush(_output_file, Z_SYNC_FLUSH);

  auto end{std::chrono::high_resolution_clock::now()};
  return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

size_t SinkGZip::FlushAndClose() {
  auto start{std::chrono::high_resolution_clock::now()};
  if (_output_file) {
    gzflush(_output_file, Z_FINISH);
    gzclose(_output_file);
    _output_file = nullptr;
  }
  auto end{std::chrono::high_resolution_clock::now()};
  return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

SinkZStd::SinkZStd(const std::filesystem::path& path, const u32 compression_level)
    : SinkBase(path), _compression_level(compression_level) {
  _zstd_context = ZSTD_createCCtx();
  if (_zstd_context == nullptr) {
    throw error::Error("Failed to create Zstandard context");
  }

  _output_file_stream.open(path, std::ios::binary);
  if (!_output_file_stream.is_open()) {
    throw error::Error("Failed to open file: {}", path.string());
  }
}

SinkZStd::~SinkZStd() { FlushAndClose(); }

size_t SinkZStd::WriteData(const SinkData& data) {
  const auto start{std::chrono::high_resolution_clock::now()};

  const auto required_buffer_size = ZSTD_compressBound(data.length);
  if (required_buffer_size > _output_buffer.size()) {
    _output_buffer.resize(required_buffer_size);
  }

  size_t const result = ZSTD_compressCCtx(_zstd_context, _output_buffer.data(), _output_buffer.size(), data.p_data,
                                          data.length, ToSigned(_compression_level));
  if (ZSTD_isError(result)) {
    throw error::Error("Zstandard error: {}. Failed to compress data.", ZSTD_getErrorName(result));
  }
  _output_file_stream.write(_output_buffer.data(), static_cast<std::streamsize>(result));
  _output_file_stream.flush();
  if (!_output_file_stream.good()) {
    throw error::Error("Failed to write to file: {}", output_file_name);
  }
  auto end{std::chrono::high_resolution_clock::now()};
  return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

size_t SinkZStd::FlushAndClose() {
  auto start{std::chrono::high_resolution_clock::now()};
  if (_output_file_stream.is_open()) {
    _output_file_stream.close();
  }
  if (_zstd_context) {
    ZSTD_freeCCtx(_zstd_context);
    _zstd_context = nullptr;
  }
  auto end{std::chrono::high_resolution_clock::now()};
  return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
}

}  // namespace xoos::demux
