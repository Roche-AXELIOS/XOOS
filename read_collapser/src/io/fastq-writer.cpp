#include "io/fastq-writer.h"

#include <zlib.h>

#include <sstream>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

namespace xoos::read_collapser {

void GzipFileDeleter::operator()(gzFile_s* gz) const {
  if (gz) {
    gzclose(gz);
  }
}

void GzipWriteFastq(gzFile_s* gz, const std::string& read, const vec<u8>& qual, const std::string& read_name) {
  std::ostringstream oss;
  oss << "@" << read_name << '\n' << read << '\n' << "+" << '\n';
  for (const auto& q : qual) {
    oss << ToPlus33Ascii(q);
  }
  oss << '\n';
  const std::string& str = oss.str();

  s32 ret = gzwrite(gz, str.c_str(), str.size());
  if (ret != static_cast<s32>(str.size())) {
    s32 err_num;
    const char* error_string = gzerror(gz, &err_num);
    if (err_num == Z_ERRNO) {
      throw error::Error("Failed to write to file: {}", std::string(error_string));
    } else {
      throw error::Error("gzip error: {}. Failed to write to file.", std::string(error_string));
    }
  }
}

void GzipWriteFastq(gzFile_s* gz, const ConsensusResult& consensus_result, const std::string& read_name) {
  std::ostringstream oss;
  oss << "@" << read_name;
  if (consensus_result.forward_per_base_depth.has_value()) {
    oss << "\tad:Z:";
    for (u32 depth : consensus_result.forward_per_base_depth.value()) {
      oss << ToPlus33Ascii(depth);
    }
  }
  if (consensus_result.reverse_per_base_depth.has_value()) {
    oss << "\tbd:Z:";
    for (u32 depth : consensus_result.reverse_per_base_depth.value()) {
      oss << ToPlus33Ascii(depth);
    }
  }
  if (consensus_result.forward_per_base_majority_count.has_value()) {
    oss << "\tam:Z:";
    for (u32 count : consensus_result.forward_per_base_majority_count.value()) {
      oss << ToPlus33Ascii(count);
    }
  }
  if (consensus_result.reverse_per_base_majority_count.has_value()) {
    oss << "\tbm:Z:";
    for (u32 count : consensus_result.reverse_per_base_majority_count.value()) {
      oss << ToPlus33Ascii(count);
    }
  }
  oss << '\n' << consensus_result.sequence << '\n' << "+" << '\n';
  for (const auto& q : consensus_result.quality_scores) {
    oss << ToPlus33Ascii(q);
  }
  oss << '\n';
  const std::string& str = oss.str();

  s32 ret = gzwrite(gz, str.c_str(), str.size());
  if (ret != static_cast<s32>(str.size())) {
    s32 err_num;
    const char* error_string = gzerror(gz, &err_num);
    if (err_num == Z_ERRNO) {
      throw error::Error("Failed to write to file: {}", std::string(error_string));
    } else {
      throw error::Error("gzip error: {}. Failed to write to file.", std::string(error_string));
    }
  }
}

GzipFilePtr OpenGzipFile(const fs::path& path, s32 compression_level) {
  create_directories(path.parent_path());
  if (compression_level > 9) {
    throw std::runtime_error("Invalid compression level: " + std::to_string(compression_level) + ", must be [0, 9]");
  }
  std::string write_mode = compression_level == 0 ? "wT" : "w" + std::to_string(compression_level);
  GzipFilePtr gz(gzopen(path.c_str(), write_mode.c_str()));
  if (!gz) {
    throw error::Error("Failed to open gzip file: '{}'", path.string());
  }
  return gz;
}

}  // namespace xoos::read_collapser
