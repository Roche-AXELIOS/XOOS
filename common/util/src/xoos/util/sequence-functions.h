#pragma once

#include <string>

namespace xoos::sequence {

// simple reverse-complement used for processing UMIs, which may be
// set to just "*" for partial reads
std::string ReverseComplement(const std::string& str);

}  // namespace xoos::sequence
