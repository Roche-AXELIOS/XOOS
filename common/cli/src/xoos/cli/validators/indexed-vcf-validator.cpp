#include "indexed-vcf-validator.h"

#include <filesystem>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <xoos/error/error.h>
#include <xoos/io/vcf/vcf-header.h>
#include <xoos/util/adopt.h>

#include "validator-util.h"
#include "xoos/io/htslib-util/htslib-util.h"

namespace xoos::cli {

void ValidateVcfAndIndex(const fs::path& path) {
  const auto vcf_file = AdoptShared(io::HtsOpen(path.c_str(), "r"));
  if (vcf_file->format.format != vcf) {
    throw error::Error("Input file is not in VCF format: '{}'", path);
  }
  // reject files that are not bgzip-compressed or are missing the EOF marker (possibly truncated)
  if (!io::HtsHasBgzfEof(vcf_file.get())) {
    throw error::Error("VCF file may be truncated: '{}'", path);
  }
  io::VcfHeader::Read(vcf_file);
  // find and open the corresponding tabix index file
  try {
    io::HtsIdxLoad(path, HTS_FMT_TBI);
  } catch (const std::runtime_error&) {
    throw error::Error("Failed to open index for VCF file: '{}'", path);
  }
}

IndexedVcfFileValidator::IndexedVcfFileValidator() {
  name_ = "INDEXED-VCF-FILE";
  func_ = CreateValidatorFunction(ValidateVcfAndIndex);
}

}  // namespace xoos::cli
