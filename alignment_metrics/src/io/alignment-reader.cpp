#include "io/alignment-reader.h"

#include <filesystem>
#include <utility>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>

namespace xoos::alignment_metrics {

AlignmentReader OpenAlignmentFile(const fs::path& location) {
  auto bam = io::HtsOpen(location, "rb");
  auto idx = io::SamIndexLoad(bam.get(), location);
  auto header = io::SamHdrRead(bam.get());
  return AlignmentReader{std::move(bam), std::move(header), std::move(idx)};
}

}  // namespace xoos::alignment_metrics
