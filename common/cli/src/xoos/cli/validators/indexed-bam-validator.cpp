#include "indexed-bam-validator.h"

#include <htslib/bgzf.h>
#include <htslib/tbx.h>

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
  {
    BGZF* bgzf = hts_get_bgzfp(bam_file.get());
    if (bgzf == nullptr) {
      throw error::Error("BAM file is not bgzip-compressed: {}", path);
    }
    if (bgzf_check_EOF(bgzf) == 0) {
      // if the EOF marker is not found, the file may be truncated
      throw error::Error("BAM file may be truncated: {}", path);
    }
  }
  const io::SamHdrPtr hdr = io::SamHdrRead(bam_file.get());
  // try reading the first record to ensure the BAM file is not empty or corrupted
  const auto record = io::Bam1Ptr(bam_init1());
  if (record == nullptr) {
    throw error::Error("Failed to allocate memory for BAM record");
  }
  const int ret = sam_read1(bam_file.get(), hdr.get(), record.get());
  if (ret == -1) {
    // EOF reached when reading the first record
    throw error::Error("BAM file is empty: '{}'", path);
  }
  if (ret < -1) {
    // cannot read the BAM file
    throw error::Error("Error reading BAM file: '{}'", path);
  }
  // find and open the corresponding index file
  io::SamIndexLoad(bam_file.get(), path.c_str());
}

IndexedBamFileValidator::IndexedBamFileValidator() {
  name_ = "INDEXED-BAM-FILE";
  func_ = CreateValidatorFunction(ValidateBamAndIndex);
}

}  // namespace xoos::cli
