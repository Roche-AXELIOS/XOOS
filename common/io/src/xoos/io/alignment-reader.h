#pragma once

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::io {

// A convenience struct to hold the bam file, header and index at once for the same file.
struct AlignmentReader {
  io::HtsFilePtr bam_file;
  io::SamHdrPtr hdr;
  io::HtsIdxPtr idx;
};

/**
 * Open a bam file in read-only mode and returns an AlignmentReader struct and update
 * the pointer to the bam file, header and index file.
 */
AlignmentReader OpenAlignmentReader(const fs::path& location);

vec<AlignmentReader> OpenAlignmentReaders(const fs::path& bam_filename, u32 reader_count);

}  // namespace xoos::io
