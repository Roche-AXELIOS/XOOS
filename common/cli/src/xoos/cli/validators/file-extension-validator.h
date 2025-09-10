#pragma once

#include <string>
#include <vector>

#include <CLI/CLI.hpp>

namespace xoos::cli {

/// Validates that the file extension is one of the allowed extensions.
class FileExtensionValidator : public CLI::Validator {
 public:
  explicit FileExtensionValidator(const std::vector<std::string>& allowed_file_extensions);
};

}  // namespace xoos::cli
