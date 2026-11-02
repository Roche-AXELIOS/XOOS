#include "bgzf-utils.h"

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>

#include <xoos/error/error.h>
#include <xoos/types/vec.h>

namespace xoos::svc {

namespace {

/// BGZF EOF marker — 28 bytes that terminate every valid .gz / .bgzf file.
/// Stored as char so it can be passed directly to ostream::write without a cast.
constexpr std::array<char, 28> kBgzfEofBlock = {
    '\x1f', '\x8b', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xff', '\x06', '\x00', '\x42', '\x43',
    '\x02', '\x00', '\x1b', '\x00', '\x03', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'};
constexpr size_t kBgzfEofSize = kBgzfEofBlock.size();

/// Append raw bytes from `src_path` to `dst`, skipping the trailing BGZF EOF marker.
/// If the file is too small to contain data beyond the EOF marker, nothing is written.
void AppendBgzfPayload(std::ofstream& dst, const fs::path& src_path) {
  const auto file_size = std::filesystem::file_size(src_path);
  if (file_size <= kBgzfEofSize) {
    return;
  }
  const auto payload_size = file_size - kBgzfEofSize;

  std::ifstream src(src_path, std::ios::binary);
  if (!src) {
    throw error::Error("Failed to open file for BGZF concatenation: {}", src_path.string());
  }

  constexpr auto kBufSize = static_cast<const size_t>(64 * 1024);
  std::array<char, kBufSize> buf{};
  size_t remaining = payload_size;
  while (remaining > 0) {
    const auto to_read = std::min(remaining, kBufSize);
    src.read(buf.data(), static_cast<std::streamsize>(to_read));
    const auto bytes_read = static_cast<size_t>(src.gcount());
    if (bytes_read == 0) {
      throw error::Error("Read failed during BGZF concatenation: {}", src_path.string());
    }
    dst.write(buf.data(), static_cast<std::streamsize>(bytes_read));
    if (!dst) {
      throw error::Error("Write failed during BGZF concatenation");
    }
    remaining -= bytes_read;
  }
}

}  // namespace

void ConcatenateBgzfFiles(const fs::path& output_path, const fs::path& header_file, const vec<fs::path>& region_files) {
  std::ofstream out(output_path, std::ios::binary);
  if (!out) {
    throw error::Error("Failed to open output file: {}", output_path.string());
  }

  // Copy header BGZF blocks (without EOF marker).
  AppendBgzfPayload(out, header_file);

  // Append each per-region data file (without EOF markers).
  for (const auto& region_file : region_files) {
    if (!std::filesystem::exists(region_file)) {
      continue;
    }
    AppendBgzfPayload(out, region_file);
  }

  // Write the final BGZF EOF marker.
  out.write(kBgzfEofBlock.data(), static_cast<std::streamsize>(kBgzfEofSize));
  if (!out) {
    throw error::Error("Failed to write BGZF EOF marker");
  }
}

}  // namespace xoos::svc
