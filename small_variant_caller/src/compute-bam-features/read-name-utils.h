#pragma once

#include <string>

#include <xoos/types/int.h>

namespace xoos::svc {

/**
 * Synopsis:
 * This file contains utility functions for parsing read names.
 */

void ParseReadName(const std::string& read_name, u32& plus_counts, u32& minus_counts, u32& family_size);

}  // namespace xoos::svc
