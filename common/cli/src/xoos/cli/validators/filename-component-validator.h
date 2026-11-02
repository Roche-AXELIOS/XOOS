#pragma once

#include <CLI/CLI.hpp>

namespace xoos::cli {

/// Validates that a value is a safe single filename component: not empty, not an absolute path, and
/// free of path separators or parent-traversal (`..`) components. Use this for options such as an
/// output prefix that are joined with a directory to build output paths, to prevent writing outside
/// the intended directory.
class FilenameComponentValidator : public CLI::Validator {
 public:
  FilenameComponentValidator();
};

}  // namespace xoos::cli
