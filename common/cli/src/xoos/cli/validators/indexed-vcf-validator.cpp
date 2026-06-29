#include "indexed-vcf-validator.h"

#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <xoos/error/error.h>
#include <xoos/io/vcf/vcf-record.h>

#include "validator-util.h"
#include "xoos/io/htslib-util/htslib-util.h"

namespace xoos::cli {

void ValidateVcfAndIndex(const fs::path& path) {
  const io::HtsFilePtr vcf_file = io::HtsOpen(path.c_str(), "r");
  if (vcf_file->format.format != vcf) {
    throw error::Error("Input file is not in VCF format: '{}'", path);
  }
  {
    BGZF* bgzf = hts_get_bgzfp(vcf_file.get());
    if (bgzf == nullptr) {
      throw error::Error("VCF file is not bgzip-compressed: {}", path);
    }
    if (bgzf_check_EOF(bgzf) == 0) {
      throw error::Error("VCF file may be truncated: {}", path);
    }
  }
  const auto hdr = io::BcfHeaderPtr(bcf_hdr_read(vcf_file.get()), bcf_hdr_destroy);
  if (hdr == nullptr) {
    throw error::Error("Failed to read VCF header: {}", path);
  }
  const auto record = io::BcfRecordPtr{bcf_init(), bcf_destroy};
  if (record == nullptr) {
    throw error::Error("Failed to allocate memory for VCF record");
  }
  const int ret = bcf_read(vcf_file.get(), hdr.get(), record.get());
  if (ret < 0) {
    throw error::Error("VCF file is empty: '{}'", path);
  }
  const io::HtsIdxPtr idx(hts_idx_load(path.c_str(), HTS_FMT_TBI));
  if (idx == nullptr) {
    throw error::Error("Failed to open index for VCF file: '{}'", path);
  }
}

IndexedVcfFileValidator::IndexedVcfFileValidator() {
  name_ = "INDEXED-VCF-FILE";
  func_ = CreateValidatorFunction(ValidateVcfAndIndex);
}

}  // namespace xoos::cli
