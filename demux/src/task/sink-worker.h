#pragma once
#include <zlib.h>
#include <zstd.h>

#include <atomic>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "io/read-record.h"

namespace xoos::demux {

constexpr u32 kMinCompressionLevel{1};
constexpr u32 kMaxGzipCompressionLevel{9};
constexpr u32 kDefaultGZipCompressionLevel{1};
constexpr u32 kMaxZstdCompressionLevel{19};
constexpr u32 kDefaultZstdCompressionLevel{1};

class FlowManager;
enum class SinkCompressionType { kNone, kGzip, kZstd };

/// @brief Sink base class for writing data to a file.
class SinkBase {
 public:
  explicit SinkBase(const std::filesystem::path& path) : output_file_name(path) {}

  virtual ~SinkBase() = default;

  // Returns the # of nanoseconds required to write out data.
  virtual size_t WriteData(const SinkData& data) = 0;
  virtual size_t FlushAndClose() = 0;

 protected:
  std::string output_file_name;
};

/// @brief Sink class for writing data to a file without applying compression.
class SinkRaw : public SinkBase {
 public:
  explicit SinkRaw(const std::filesystem::path& path);
  ~SinkRaw() override;

  size_t WriteData(const SinkData& data) override;
  size_t FlushAndClose() override;

 private:
  std::ofstream _output_file_stream;
};

/// @brief Sink class for writing data to a file with GZip compression.
class SinkGZip : public SinkBase {
 public:
  SinkGZip(const std::filesystem::path& path, size_t compression_level);
  ~SinkGZip() override;

  size_t WriteData(const SinkData& data) override;
  size_t FlushAndClose() override;

 private:
  gzFile _output_file{nullptr};
  std::atomic<int> _nr_writers{0};
  std::mutex _write_mutex;
};

/// @brief Sink class for writing data to a file with ZStandard compression.
class SinkZStd : public SinkBase {
 public:
  SinkZStd(const std::filesystem::path& path, u32 compression_level);
  ~SinkZStd() override;

  size_t WriteData(const SinkData& data) override;
  size_t FlushAndClose() override;

 private:
  ZSTD_CCtx* _zstd_context{nullptr};
  std::ofstream _output_file_stream;
  std::vector<char> _output_buffer;
  u32 _compression_level{1};
};

/// @brief Sink class for writing data to a Null sink (i.e. no output).
class SinkNull : public SinkBase {
 public:
  explicit SinkNull(const std::filesystem::path& path) : SinkBase(path) {}

  ~SinkNull() override = default;

  size_t WriteData(const SinkData&) override { return 5ul; }  // Return 5 ns as time

  size_t FlushAndClose() override { return 5ul; }
};
}  // namespace xoos::demux
