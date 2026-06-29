#include "bed-validator.h"

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-ptr.h>

#include "validator-util.h"
#include "xoos/io/htslib-util/htslib-util.h"

namespace xoos::cli {

void ValidateBed(const fs::path& path) {
  const io::HtsFilePtr hts_file = io::HtsOpen(path.c_str(), "r");
  if (hts_file->format.format != bed) {
    throw error::Error("Input file is not in BED format: '{}'", path);
  }
}

BedFileValidator::BedFileValidator() {
  name_ = "BED-FILE";
  func_ = CreateValidatorFunction(ValidateBed);
}

}  // namespace xoos::cli
