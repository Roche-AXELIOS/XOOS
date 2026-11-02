#include "bed-validator.h"

#include <xoos/error/error.h>
#include <xoos/io/fast-bed-parser.h>

#include "validator-util.h"

namespace xoos::cli {

void ValidateBed(const fs::path& path) {
  bool found_region = false;
  io::ForEachBedRegionInFileFast(path, [&found_region](io::BedRegion&) {
    found_region = true;
    return false;
  });
  if (!found_region) {
    throw error::Error("BED file is empty or does not contain any valid regions");
  }
}

BedFileValidator::BedFileValidator() {
  name_ = "BED-FILE";
  func_ = CreateValidatorFunction(ValidateBed);
}

}  // namespace xoos::cli
