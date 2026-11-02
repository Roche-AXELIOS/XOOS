#pragma once
#include <xoos/types/fs.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {

using SeqLenMap = StrUnorderedMap<size_t>;
SeqLenMap ParseFai(const fs::path& fname);
StrUnorderedMap<size_t> GetContigOrder(const fs::path& fname);
vec<std::string> GetContigsInOrder(const fs::path& fname);
}  // namespace xoos::cnc
