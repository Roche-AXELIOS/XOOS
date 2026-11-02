#include "indexed-bam-validator.h"

#include <filesystem>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include "validator-util.h"
#include "xoos/error/error.h"
#include "xoos/io/htslib-util/htslib-util.h"

namespace xoos::cli {

void ValidateBamAndIndex(const fs::path& path) {
  // open the BAM file with htslib
  const io::HtsFilePtr bam_file = io::HtsOpen(path.c_str(), "r");
  if (bam_file->format.format != bam) {
    throw error::Error("Input file is not in BAM format: '{}'", path);
  }
  // reject files that are not bgzip-compressed or are missing the EOF marker (possibly truncated)
  if (!io::HtsHasBgzfEof(bam_file.get())) {
    throw error::Error("BAM file may be truncated: '{}'", path);
  }
  io::SamHdrRead(bam_file.get());
  // find and open the corresponding index file
  io::SamIndexLoad(bam_file.get(), path.c_str());
}

IndexedBamFileValidator::IndexedBamFileValidator() {
  name_ = "INDEXED-BAM-FILE";
  func_ = CreateValidatorFunction(ValidateBamAndIndex);
}

}  // namespace xoos::cli
