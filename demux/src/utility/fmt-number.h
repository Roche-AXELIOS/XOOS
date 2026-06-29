#pragma once

#include <string>

namespace xoos::demux {

/**
 * Given a number, return a string formatted with the number part <1000
 * and an appropriate metric suffix:
 *   kilo	k	1000
 *   mega	M	1000000
 *   giga	G	1000000000
 *   tera	T	1000000000000
 *   peta	P	1000000000000000
 *   exa	E	1000000000000000000
 */
std::string FormatWithMetricSuffix(double value);

}  // namespace xoos::demux
