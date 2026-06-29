#pragma once

#include <optional>
#include <string>

#include "sequence/loci-range.h"

namespace xoos::demux {
std::string_view Substr(const std::string_view& str, std::string::size_type start_pos, std::string::size_type end_pos);

std::string Reverse(const std::string_view& str);

std::optional<LociRange> FindExactMatch5p(const char* search_seq, uint search_seq_length, int max_wiggle_right,
                                          uint start_pos, const char* read_seq, uint seq_length);

std::optional<LociRange> FindExactMatch3p(const char* search_seq, uint search_seq_len, int max_wiggle_left,
                                          uint start_pos, const char* read_seq, uint seq_length);
}  // namespace xoos::demux
