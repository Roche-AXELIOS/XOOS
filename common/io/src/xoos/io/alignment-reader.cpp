#include "xoos/io/alignment-reader.h"

#include <filesystem>
#include <string>
#include <utility>

#include <htslib/hts.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::io {

AlignmentReader OpenAlignmentReader(const fs::path& location) {
  auto bam_file = io::HtsOpen(location, "rb");
  if (bam_file->format.format != htsExactFormat::sam && bam_file->format.format != htsExactFormat::bam) {
    throw error::Error(
        "Unsupported alignment file format when opening '{}': detected HTS format code {} (supported formats: SAM, "
        "BAM)",
        location.string(),
        std::to_string(static_cast<s32>(bam_file->format.format)));
  }

  auto idx = io::HtsIdxLoad(location, HTS_FMT_BAI);
  auto hdr = io::SamHdrRead(bam_file.get());
  return AlignmentReader{.bam_file = std::move(bam_file), .hdr = std::move(hdr), .idx = std::move(idx)};
}

vec<AlignmentReader> OpenAlignmentReaders(const fs::path& bam_filename, const u32 reader_count) {
  vec<AlignmentReader> alignment_readers;
  alignment_readers.reserve(reader_count);
  for (u32 i = 0; i < reader_count; ++i) {
    alignment_readers.emplace_back(OpenAlignmentReader(bam_filename));
  }
  return alignment_readers;
}

}  // namespace xoos::io
