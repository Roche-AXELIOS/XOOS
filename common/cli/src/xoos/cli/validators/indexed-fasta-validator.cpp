#include "indexed-fasta-validator.h"

#include <filesystem>

#include <htslib/faidx.h>
#include <htslib/hts.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-ptr.h>

#include "validator-util.h"
#include "xoos/io/htslib-util/htslib-util.h"

namespace xoos::cli {

void ValidateFastaAndIndex(const fs::path& path) {
  const io::HtsFilePtr hts_file = io::HtsOpen(path.c_str(), "r");
  if (hts_file->format.format != fasta_format) {
    throw error::Error("Input file is not in FASTA format: '{}'", path);
  }
  // Validator only succeeds if the index already exists (FaiLoad disables automatic index creation).
  io::FaiLoad(path);
}

IndexedFastaFileValidator::IndexedFastaFileValidator() {
  name_ = "INDEXED-FASTA-FILE";
  func_ = CreateValidatorFunction(ValidateFastaAndIndex);
}

}  // namespace xoos::cli
