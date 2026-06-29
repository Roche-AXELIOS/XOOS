#pragma once

#include <zlib.h>

#include <concepts>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "consensus/consensus-result.h"

namespace xoos::read_collapser {

// Cap at 93 and add 33 to fit in ASCII printable range
template <std::unsigned_integral T>
static inline char ToPlus33Ascii(const T val) {
  return static_cast<char>(std::min(val, static_cast<T>(93)) + static_cast<T>(33));
}

struct GzipFileDeleter {
  void operator()(gzFile_s* gz) const;
};

using GzipFilePtr = std::unique_ptr<gzFile_s, GzipFileDeleter>;

void GzipWriteFastq(gzFile_s* gz, const std::string& read, const vec<u8>& qual, const std::string& read_name);

void GzipWriteFastq(gzFile_s* gz, const ConsensusResult& consensus_result, const std::string& read_name);

GzipFilePtr OpenGzipFile(const fs::path& path, s32 compression_level = 0);

}  // namespace xoos::read_collapser
