#pragma once

#include <algorithm>
#include <optional>
#include <string>

namespace xoos::util::container {

template <typename C, typename K = C::key_type, typename V = C::value_type::second_type>
std::optional<V> Find(const C& elements, const K& k) {
  auto it = elements.find(k);
  if (it == std::end(elements)) {
    return std::nullopt;
  }
  return it->second;
}

template <typename C, typename T>
bool Contains(const C& elements, const T& value) {
  return std::find(std::cbegin(elements), std::cend(elements), value) != std::cend(elements);
}

}  // namespace xoos::util::container
