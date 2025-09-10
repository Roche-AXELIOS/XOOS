#include "gzip.h"

#include <fstream>

#include <xoos/error/error.h>
#include <xoos/types/vec.h>

#include "util/compress/ios-util.h"

namespace xoos::svc {

void CompressGzip(const fs::path& input_path, const fs::path& output_path, size_t buffer_size, int compression_level) {
  auto in = OpenInput(input_path, std::ios::binary);

  auto out = GzOpen(output_path, "wb");
  GzSetParams(out.get(), compression_level, Z_DEFAULT_STRATEGY);

  vec<char> buffer(buffer_size);
  while (in) {
    const auto bytes_read = Read(in, buffer);
    GzWrite(out.get(), buffer.data(), bytes_read);
  }

  RequireIoNotBad(input_path, in);
}

void CompressGzip(const std::string& input_content, const fs::path& output_path, int compression_level) {
  auto out = GzOpen(output_path, "wb");
  GzSetParams(out.get(), compression_level, Z_DEFAULT_STRATEGY);
  GzWrite(out.get(), input_content.data(), input_content.size());
}

void DecompressGzip(const fs::path& input_path, const fs::path& output_path, size_t buffer_size) {
  auto in = GzOpen(input_path, "rb");

  auto out = OpenOutput(output_path, std::ios::binary);

  vec<char> buffer(buffer_size);
  while (true) {
    const auto bytes_read = GzRead(in.get(), buffer.data(), buffer.size());
    if (bytes_read > 0) {
      out.write(buffer.data(), bytes_read);
    }
    if (std::cmp_less(bytes_read, buffer.size())) {
      break;
    }
  }

  RequireIoNotBad(output_path, out);
}

std::string DecompressGzip(const fs::path& input_path, size_t buffer_size) {
  auto in = GzOpen(input_path, "rb");

  std::string out;
  const auto estimated_compression_ratio = 3;
  out.reserve(fs::file_size(input_path) * estimated_compression_ratio);

  vec<char> buffer(buffer_size);
  while (true) {
    const auto bytes_read = GzRead(in.get(), buffer.data(), buffer.size());
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

void GzWrite(gzFile file, const char* data, size_t size) {
  if (gzwrite(file, data, size) == 0) {
    GzError(file);
  }
}

void GzSetParams(gzFile file, int level, int strategy) {
  if (gzsetparams(file, level, strategy) == -1) {
    GzError(file);
  }
}

int GzRead(gzFile file, char* buffer, size_t size) {
  int bytes_read = gzread(file, buffer, size);
  if (bytes_read <= 0) {
    GzError(file);
  }
  return bytes_read;
}

}  // namespace xoos::svc
