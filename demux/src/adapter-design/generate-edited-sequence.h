#pragma once

#include <xoos/types/int.h>

#include <functional>
#include <string>

namespace xoos::demux {
/**
 * A callback function that will be invoked for each edited sequence,
 * the first argument is the edit distance of the provided sequence, the second is
 * the sequence itself
 */
using EditedSequenceCb = std::function<void(u32, const std::string&)>;

/**
 * Generate all edited sequence derived from a given sequence up to a maximum edit distance.
 * @param alphabet all possible characters that can be used in an edited sequence
 * @param sequence the original sequence to edit
 * @param max_edit_distance the largest edit distance for which to produce edited sequences
 * @param cb all callback function invoked with each edited sequence
 */
void GenerateEditedSequence(const std::string_view& alphabet, const std::string_view& sequence, u32 max_edit_distance,
                            const EditedSequenceCb& cb);
}  // namespace xoos::demux
