#pragma once

#include <gsl/gsl>  // NOLINT

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/fs.h>

namespace xoos::alignment_metrics {

// A convenience struct to hold the bam file, header and index at once for the same file.
struct AlignmentReader {
  io::HtsFilePtr bam;
  io::SamHdrPtr header;
  io::HtsIdxPtr idx;
};

/**
 * Open a bam file in read-only mode and returns an AlignmentReader struct and update
 * the pointer to the bam file, header and index file.
 */
AlignmentReader OpenAlignmentFile(const fs::path& location);

}  // namespace xoos::alignment_metrics
