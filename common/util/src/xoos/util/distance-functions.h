#pragma once

#include <string>

namespace xoos::distance {

/**
 * Calculate the hamming distance between two strings, exit early once a distance threshold is
 * exceeded.
 */
bool HammingDistanceThreshold(const std::string& source, const std::string& target, size_t distance);

/**
 * Calculate the levenshtein distance between two strings, exit early once a distance threshold is
 * exceeded.
 */
bool LevenshteinDistanceThreshold(const std::string& source, const std::string& target, size_t distance);

/**
 * First calculate the hamming distance between source and target, if the distance threshold is exceeded then
 * exit early with false, if not calculate the levenshtein distance, if the distance threshold is exceeded then exit
 * early with false.
 *
 * Has special behaviour to handle UMI strings, in particular if either string is '*' exit early with false.
 */
bool IsUmiDistanceLeq(const std::string& source, const std::string& target, size_t distance);

}  // namespace xoos::distance
