#pragma once

#include <memory>
#include <vector>

namespace xoos {

template <typename Tp, typename Alloc = std::allocator<Tp>>
using vec = std::vector<Tp, Alloc>;

}  // namespace xoos
