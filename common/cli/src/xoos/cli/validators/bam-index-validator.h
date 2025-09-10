#pragma once

#include <CLI/Validators.hpp>

namespace xoos::cli {

// A validator to check that a BAM index exists for a given BAM file.
struct BamIndexValidator : CLI::Validator {
  BamIndexValidator();
};

static const BamIndexValidator kBamIndex;

}  // namespace xoos::cli
