#pragma once

#include <vector>

namespace xoos::demux {

template <typename T>
std::vector<T> Range(T start, T end, T step) {
  auto result = std::vector<T>{};
  for (auto i = start; i < end; i += step) {
    result.emplace_back(i);
  }
  return result;
}

}  // namespace xoos::demux
