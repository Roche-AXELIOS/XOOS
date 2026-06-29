#include "indexed-fasta-validator.h"

#include <xoos/error/error.h>

#include "validator-util.h"
#include "xoos/io/htslib-util/htslib-util.h"

namespace xoos::cli {

void ValidateFastaAndIndex(const fs::path& path) {
  const io::HtsFilePtr hts_file = io::HtsOpen(path.c_str(), "r");
  if (hts_file->format.format != fasta_format) {
    throw error::Error("Input file is not in FASTA format: '{}'", path);
  }
  // Use fai_load3 with creation flag 0 to prevent automatic index creation; validator only succeeds if index exists.
  const io::FaIdxPtr fai(fai_load3(path.c_str(), nullptr, nullptr, 0));
  if (fai == nullptr) {
    throw error::Error("Failed to open index for FASTA file: '{}'", path);
  }
}

IndexedFastaFileValidator::IndexedFastaFileValidator() {
  name_ = "INDEXED-FASTA-FILE";
  func_ = CreateValidatorFunction(ValidateFastaAndIndex);
}

}  // namespace xoos::cli
