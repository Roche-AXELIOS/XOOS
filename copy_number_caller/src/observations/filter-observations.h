#pragma once
#include <string>

#include <xoos/types/vec.h>

#include "observations/observations.h"

namespace xoos::cnc {
Observations GetOverlappingObservations(const Observations& left, const Observations& right);
Observations GetOverlappingObservations(const vec<std::string>& left, const Observations& right);
}  // namespace xoos::cnc
