#include "compress.h"

#include <xoos/compress/gzip.h>
#include <xoos/compress/zstandard.h>
#include <xoos/error/error.h>

namespace xoos::compress {

static const std::string kGzipExtension = ".gz";
constexpr size_t kGzipDefaultBufferSize = 4096;
constexpr int kGzipDefaultCompressionLevel = 6;

static const std::string kZstdExtension = ".zst";
constexpr int kZstdDefaultCompressionLevel = 3;
constexpr int kZstdDefaultThreads = 1;

bool IsCompressed(const fs::path& path) {
  return path.extension() == kGzipExtension || path.extension() == kZstdExtension;
}

void Compress(const fs::path& input_path, const fs::path& output_path) {
  if (output_path.extension() == kGzipExtension) {
    compress::CompressGzip(input_path, output_path, kGzipDefaultBufferSize, kGzipDefaultCompressionLevel);
  } else if (output_path.extension() == kZstdExtension) {
    compress::CompressZstd(input_path, output_path, kZstdDefaultCompressionLevel, kZstdDefaultThreads);
  } else {
    throw error::Error("Unsupported compression format for file: {}", output_path);
  }
}

void Compress(const std::string& input_content, const fs::path& output_path) {
  if (output_path.extension() == kGzipExtension) {
    compress::CompressGzip(input_content, output_path, kGzipDefaultCompressionLevel);
  } else if (output_path.extension() == kZstdExtension) {
    compress::CompressZstd(input_content, output_path, kZstdDefaultCompressionLevel);
  } else {
    throw error::Error("Unsupported compression format for file: {}", output_path);
  }
}

void Decompress(const fs::path& input_path, const fs::path& output_path) {
  if (input_path.extension() == kGzipExtension) {
    compress::DecompressGzip(input_path, output_path, kGzipDefaultBufferSize);
  } else if (input_path.extension() == kZstdExtension) {
    compress::DecompressZstd(input_path, output_path);
  } else {
    throw error::Error("Unsupported compression format for file: {}", input_path);
  }
}

std::string Decompress(const fs::path& input_path) {
  if (input_path.extension() == kGzipExtension) {
    return compress::DecompressGzip(input_path, kGzipDefaultBufferSize);
  } else if (input_path.extension() == kZstdExtension) {
    return compress::DecompressZstd(input_path);
  } else {
    throw error::Error("Unsupported compression format for file: {}", input_path);
  }
}

}  // namespace xoos::compress
