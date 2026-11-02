#include "xoos/cli/validators/filename-component-validator.h"

#include <string>

#include <fmt/format.h>

#include <xoos/types/fs.h>

#include "CLI/CLI.hpp"

namespace xoos::cli {

FilenameComponentValidator::FilenameComponentValidator() : Validator("FILENAME_COMPONENT") {
  func_ = [](std::string& value) -> std::string {
    if (value.empty()) {
      return "File name must not be empty.";
    }
    // Reject any separator from either platform so behavior is consistent regardless of host OS.
    if (value.contains('/') || value.contains('\\')) {
      return fmt::format("File name must be a single filename component and must not contain path separators: {}",
                         value);
    }
    const fs::path path{value};
    if (path.is_absolute()) {
      return fmt::format("File name must be a relative filename component, not an absolute path: {}", value);
    }
    if (value == "." || value == "..") {
      return fmt::format("File name must not be a path-traversal component ('.' or '..'): {}", value);
    }
    return {};  // Valid
  };
}

}  // namespace xoos::cli
