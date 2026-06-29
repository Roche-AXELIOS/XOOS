#pragma once

#include <cstddef>

#include <xoos/types/float.h>

namespace xoos::read_collapser {

// Represents the indices corresponding to different bases and special symbols used in the consensus matrix.
enum class BaseIndex : size_t {
  kIndexA = 0,
  kIndexC = 1,
  kIndexG = 2,
  kIndexT = 3,
  kIndexGap = 4,
  kIndexN = 5,
  kIndexP = 6
};

// Total number of base indices including A, C, G, T, Gap, N, and P.
static constexpr size_t kBaseIndexCount = 7;

static constexpr char kBaseA = 'A';
static constexpr char kBaseC = 'C';
static constexpr char kBaseG = 'G';
static constexpr char kBaseT = 'T';
/**
 * 'N' represents an ambiguous base in the consensus matrix.
 *
 * It does not count towards the total depth of a position and does not
 * participate in majority voting.
 */
static constexpr char kBaseN = 'N';
/**
 * 'P' is used to indicate a partial base in the consensus matrix
 * It is used to indicate that the current read does not cover the position
 * but some other read does. For example,
 *
 * Read 1: ACGTACGT
 * Read 2: ACGTAC
 *
 * Consensus Matrix:
 * ACGTACGT
 * ACGTACPP
 *
 * 'P' does not count towards the total depth of a position and does not
 * participate in majority voting.
 */
static constexpr char kBaseP = 'P';
/**
 * A gap is used to indicate that a position is deleted base in the current read
 * or that some other read has an insertion at this position. For example,
 *
 * Read 1: ACGTACGT
 * Read 2: ACGTCGT
 *
 * Consensus Matrix:
 * ACGTACGT
 * ACGT-CGT
 *
 * A gap counts towards the total depth of a position and participates in majority voting.
 */
static constexpr char kBaseGap = '-';

/**
 * Converts a base character to its corresponding BaseIndex.
 *
 * @param base The base character to convert. Expected values are:
 *             'A', 'C', 'G', 'T', 'N', 'P', or '-'.
 * @return The corresponding BaseIndex value for the given base character.
 * Note that invalid base characters will return kIndexGap.
 */
BaseIndex ToIndex(char base);

/**
 * Converts a base character to its corresponding size_t index.
 *
 * @param base The base character to convert. Expected values are:
 *             'A', 'C', 'G', 'T', 'N', 'P', or '-'.
 * @return The size_t index corresponding to the given base character.
 *         The index values are consistent with the BaseIndex.
 * Note that invalid base characters will return the size_t index for kIndexGap.
 */
size_t ToSizeTIndex(char base);

/**
 * Converts a BaseIndex value to its corresponding base character.
 *
 * @param index The BaseIndex value to convert.
 * @return The base character corresponding to the given BaseIndex value.
 *         Returns one of 'A', 'C', 'G', 'T', 'N', 'P', or '-'.
 * Note that invalid indices will return kBaseGap.
 */
char ToBase(BaseIndex index);

}  // namespace xoos::read_collapser
