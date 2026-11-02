#include "task/sink-worker.h"

#include <array>
#include <limits>

#include <fmt/format.h>

#include <xoos/error/error.h>

#include "utility/stop-watch.h"

namespace xoos::demux {

SinkRaw::SinkRaw(const std::filesystem::path& path) : SinkBase(path) {
  _output_file_stream.open(path, std::ios::binary);
  if (!_output_file_stream.is_open()) {
    throw error::Error("Failed to open file: {}", path.string());
  }
}

SinkRaw::~SinkRaw() {
  SinkRaw::FlushAndClose();
}

size_t SinkRaw::WriteData(const SinkData& data) {
  const StopWatch sw;
  _output_file_stream.write(data.p_data, static_cast<std::streamsize>(data.length));
  if (!_output_file_stream.good()) {
    throw error::Error("Failed to write to file: {}", output_file_name);
  }
  _output_file_stream.flush();
  return ToUnsigned(sw.ElapsedTime());
}

size_t SinkRaw::FlushAndClose() {
  const auto start{std::chrono::high_resolution_clock::now()};
  if (_output_file_stream.is_open()) {
    _output_file_stream.flush();
    _output_file_stream.close();
  }
  const auto end{std::chrono::high_resolution_clock::now()};
  return ToUnsigned(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
}

SinkGZip::SinkGZip(const std::filesystem::path& path, size_t compression_level) : SinkBase(path) {
  const std::string write_mode = fmt::format("w{}", compression_level);
  _output_file = gzopen(path.c_str(), write_mode.c_str());
  if (_output_file == Z_NULL) {
    throw error::Error("Failed to open file: '{}'", path.string());
  }
}

SinkGZip::~SinkGZip() {
  SinkGZip::FlushAndClose();
}

size_t SinkGZip::WriteData(const SinkData& data) {
  const auto start{std::chrono::high_resolution_clock::now()};

  // gzwrite takes a 32-bit unsigned length; reject buffers that would not fit before writing anything, so an
  // oversized batch fails cleanly instead of being silently truncated into a corrupt output file.
  if (data.length > std::numeric_limits<u32>::max()) {
    throw error::Error("Buffer of {} bytes exceeds the maximum gzwrite length. Failed to write to file: {}",
                       data.length,
                       output_file_name);
  }

  // gzwrite returns the number of bytes written as int (0 on error).
  // NOSONAR(S5356): data.p_data (char*) is passed to zlib's const void* parameter at the C-API boundary.
  const auto bytes_written = gzwrite(_output_file, data.p_data, static_cast<u32>(data.length));  // NOSONAR(S5356)
  const size_t nr_written = ToUnsigned(bytes_written);
  if (nr_written != data.length) {
    s32 err_num = 0;
    const char* error_string = gzerror(_output_file, &err_num);
    if (err_num == Z_ERRNO) {
      throw error::Error("Failed to write to file: {}. System error: {}", output_file_name, error_string);
    }
    throw error::Error("gzip error: {}. Failed to write to file: {}", std::string(error_string), output_file_name);
  }
  gzflush(_output_file, Z_SYNC_FLUSH);

  const auto end{std::chrono::high_resolution_clock::now()};
  return ToUnsigned(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
}

size_t SinkGZip::FlushAndClose() {
  const auto start{std::chrono::high_resolution_clock::now()};
  if (_output_file != nullptr) {
    gzflush(_output_file, Z_FINISH);
    gzclose(_output_file);
    _output_file = nullptr;
  }
  const auto end{std::chrono::high_resolution_clock::now()};
  return ToUnsigned(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
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

SinkZStd::~SinkZStd() {
  SinkZStd::FlushAndClose();
}

size_t SinkZStd::WriteData(const SinkData& data) {
  const auto start{std::chrono::high_resolution_clock::now()};

  const auto required_buffer_size = ZSTD_compressBound(data.length);
  if (required_buffer_size > _output_buffer.size()) {
    _output_buffer.resize(required_buffer_size);
  }

  // NOSONAR(S5356): the buffer (char*) and data.p_data (char*) are passed to ZSTD's void*/const void*
  // parameters at the C-API boundary; there is no related-pointer alternative.
  const size_t result = ZSTD_compressCCtx(_zstd_context,
                                          _output_buffer.data(),
                                          _output_buffer.size(),
                                          data.p_data,
                                          data.length,
                                          ToSigned(_compression_level));
  if (ZSTD_isError(result) != 0U) {
    throw error::Error("Zstandard error: {}. Failed to compress data.", ZSTD_getErrorName(result));
  }
  _output_file_stream.write(_output_buffer.data(), static_cast<std::streamsize>(result));
  _output_file_stream.flush();
  if (!_output_file_stream.good()) {
    throw error::Error("Failed to write to file: {}", output_file_name);
  }
  _data_written = true;
  const auto end{std::chrono::high_resolution_clock::now()};
  return ToUnsigned(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
}

size_t SinkZStd::FlushAndClose() {
  const auto start{std::chrono::high_resolution_clock::now()};
  if (_output_file_stream.is_open()) {
    // A sink that received no data (empty-output case) would otherwise leave a zero-byte file, which is not a
    // valid Zstandard frame. Emit an empty frame so the output is well-formed and decodable by strict readers, matching
    // the gzip sink which writes a valid empty member.
    if (!_data_written) {
      std::array<char, ZSTD_COMPRESSBOUND(0)> empty_frame{};
      const size_t result = ZSTD_compressCCtx(
          _zstd_context, empty_frame.data(), empty_frame.size(), nullptr, 0, ToSigned(_compression_level));
      if (ZSTD_isError(result) == 0U) {
        _output_file_stream.write(empty_frame.data(), static_cast<std::streamsize>(result));
        _output_file_stream.flush();
      }
    }
    _output_file_stream.close();
  }
  if (_zstd_context != nullptr) {
    ZSTD_freeCCtx(_zstd_context);
    _zstd_context = nullptr;
  }
  const auto end{std::chrono::high_resolution_clock::now()};
  return ToUnsigned(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
}

}  // namespace xoos::demux
