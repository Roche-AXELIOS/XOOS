#include "xoos/cli/validators/bam-index-validator.h"

#include <filesystem>

namespace xoos::cli {

namespace fs = std::filesystem;

BamIndexValidator::BamIndexValidator() noexcept {
  name_ = "BAM INDEX(existing)";
  func_ = [](const std::string& input) {
    const auto bam_index = input + ".bai";
    if (!fs::exists(bam_index)) {
      return "Path does not exist " + bam_index;
    }
    return std::string{};
  };
}
}  // namespace xoos::cli
