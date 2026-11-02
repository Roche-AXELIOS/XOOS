#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <xoos/types/float.h>
#include <xoos/types/vec.h>

namespace xoos::gc_correct {

// Genomic region: 0-based half-open [start, end).
struct GenomicInterval {
  std::string contig;
  // 0-based inclusive start position
  size_t start;
  // 0-based exclusive end position
  size_t end;

  std::string ToString() const {
    return fmt::format("{}:{}-{}", contig, start, end);
  }
};

struct CoverageRecord {
  GenomicInterval region;
  f64 read_count = 0;
  f64 total_coverage = 0;
  f64 average_coverage = 0;
};

using CoverageRecords = std::vector<CoverageRecord>;

}  // namespace xoos::gc_correct
