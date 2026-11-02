#include "xoos/cli/validators/file-extension-validator.h"

#include <string>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <xoos/util/string-functions.h>

#include "CLI/CLI.hpp"

namespace xoos::cli {
FileExtensionValidator::FileExtensionValidator(const std::vector<std::string>& allowed_file_extensions)
    : Validator(fmt::format("FILE_EXTENSION({})", fmt::join(allowed_file_extensions, ", "))) {
  func_ = [allowed_file_extensions](std::string& filename) -> std::string {
    if (!string::EndsWithOneOf(filename, allowed_file_extensions)) {
      return fmt::format(
          "File must have one of the following extensions: {}: {}", fmt::join(allowed_file_extensions, ", "), filename);
    }
    return {};  // Valid
  };
}
}  // namespace xoos::cli
