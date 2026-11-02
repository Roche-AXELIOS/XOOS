#include "gzip.h"

#include <zlib.h>

#include <cstddef>
#include <fstream>
#include <ios>
#include <limits>
#include <utility>

#include <xoos/error/error.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/ios-util.h>

namespace xoos::compress {

void CompressGzip(const fs::path& input_path, const fs::path& output_path, size_t buffer_size, int compression_level) {
  auto in = util::OpenInput(input_path, std::ios::binary);

  auto out = GzOpen(output_path, "wb");
  GzSetParams(out.get(), compression_level, Z_DEFAULT_STRATEGY);

  vec<char> buffer(buffer_size);
  while (in) {
    const auto bytes_read = util::Read(in, buffer);
    GzWrite(out.get(), buffer.data(), bytes_read);
  }

  util::RequireIoNotBad(input_path, in);
}

void CompressGzip(const std::string& input_content, const fs::path& output_path, int compression_level) {
  auto out = GzOpen(output_path, "wb");
  GzSetParams(out.get(), compression_level, Z_DEFAULT_STRATEGY);
  GzWrite(out.get(), input_content.data(), input_content.size());
}

void DecompressGzip(const fs::path& input_path, const fs::path& output_path, size_t buffer_size) {
  auto in = GzOpen(input_path, "rb");

  auto out = util::OpenOutput(output_path, std::ios::binary);

  vec<char> buffer(buffer_size);
  while (true) {
    const auto bytes_read = GzRead(in.get(), buffer);
    if (bytes_read > 0) {
      out.write(buffer.data(), bytes_read);
    }
    if (std::cmp_less(bytes_read, buffer.size())) {
      break;
    }
  }

  util::RequireIoNotBad(output_path, out);
}

std::string DecompressGzip(const fs::path& input_path, size_t buffer_size) {
  auto in = GzOpen(input_path, "rb");

  std::string out;
  const auto estimated_compression_ratio = 3;
  out.reserve(fs::file_size(input_path) * estimated_compression_ratio);

  vec<char> buffer(buffer_size);
  while (true) {
    const auto bytes_read = GzRead(in.get(), buffer);
    if (bytes_read > 0) {
      out.append(buffer.data(), bytes_read);
    }
    if (std::cmp_less(bytes_read, buffer.size())) {
      break;
    }
  }
  return out;
}

void GzFileDeleter::operator()(gzFile file) const {
  if (file != nullptr) {
    gzclose(file);
  }
}

GzFilePtr GzOpen(const fs::path& path, const std::string& mode) {
  gzFile file = gzopen(path.c_str(), mode.c_str());
  if (file == nullptr) {
    throw error::Error("Failed to open gzip file: '{}'", path);
  }
  return GzFilePtr{file};
}

void GzError(gzFile file) {
  int err{};
  const auto* msg = gzerror(file, &err);
  if (err != Z_OK) {
    throw error::Error("Gzip error: {}", msg);
  }
}

void GzWrite(gzFile file, const char* const data, const u32 size) {
  if (gzwrite(file, data, size) == 0) {  // NOSONAR(S5356) — zlib C API requires void*
    GzError(file);
  }
}

void GzSetParams(gzFile file, int level, int strategy) {
  if (gzsetparams(file, level, strategy) == -1) {
    GzError(file);
  }
}

s32 GzRead(gzFile file, char* const buffer, const u32 size) {
  const auto bytes_read = gzread(file, buffer, size);  // NOSONAR(S5356) — zlib C API requires void*
  if (bytes_read < 0) {
    GzError(file);
  }
  return bytes_read;
}

void GzWrite(gzFile file, const char* const data, const size_t size) {
  constexpr auto kMaxSize = std::numeric_limits<u32>::max();
  if (size > kMaxSize) {
    throw error::Error("GzWrite size {} exceeds maximum of {}", size, kMaxSize);
  }
  GzWrite(file, data, static_cast<u32>(size));
}

s32 GzRead(gzFile file, char* const buffer, const size_t size) {
  constexpr auto kMaxSize = std::numeric_limits<u32>::max();
  if (size > kMaxSize) {
    throw error::Error("GzRead size {} exceeds maximum of {}", size, kMaxSize);
  }
  return GzRead(file, buffer, static_cast<u32>(size));
}

s32 GzRead(gzFile file, vec<char>& buffer) {
  return GzRead(file, buffer.data(), buffer.size());
}

}  // namespace xoos::compress
