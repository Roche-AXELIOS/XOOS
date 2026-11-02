#pragma once

#include <string>

#include <xoos/types/fs.h>

namespace xoos::cnc {

// this excludes: unmapped, secondary, QC fail and duplicate reads
const std::string kCalculateCoverageDefaultExcludeFlags = "1796";
const bool kCalculateCoverageDefaultIgnoreDN = false;

struct CalculateCoverageOptions {
  std::string exclude_flags = kCalculateCoverageDefaultExcludeFlags;
  bool ignore_DN = kCalculateCoverageDefaultIgnoreDN;
};

}  // namespace xoos::cnc
