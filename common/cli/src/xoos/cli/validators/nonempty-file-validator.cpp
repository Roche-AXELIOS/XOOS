#include "nonempty-file-validator.h"

#include <xoos/error/error.h>

#include "validator-util.h"

namespace xoos::cli {

void ValidateNonEmptyFile(const fs::path& path) {
  if (fs::file_size(path) == 0) {
    throw error::Error("File is empty: {}", path);
  }
}

NonEmptyFileValidator::NonEmptyFileValidator() {
  name_ = "NONEMPTY-FILE";
  func_ = CreateValidatorFunction(ValidateNonEmptyFile);
}

}  // namespace xoos::cli
