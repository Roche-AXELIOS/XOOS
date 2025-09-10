#pragma once

#include <cstdint>
#include <functional>

namespace xoos::util::hash {

// Below hash functions are very similar to:
// D. Eastlake, “The FNV Non-Cryptographic Hash Algorithm,” Internet-Draft draft-eastlake-fnv-31,
// Internet Engineering Task Force, 2024. [Online]. Available:
// https://www.ietf.org/archive/id/draft-eastlake-fnv-31.html.

template <typename... Ts>
std::uint64_t HashCombine(Ts... values) {
  constexpr std::uint64_t kPrime{0x100000001b3};
  std::uint64_t result{0xcbf29ce484222325};
  for (const auto& value : {values...}) {
    result = (result * kPrime) ^ value;
  }
  return result;
}

template <typename... Ts>
std::uint64_t Hash(Ts... values) {
  constexpr std::uint64_t kPrime{0x100000001b3};
  std::uint64_t result{0xcbf29ce484222325};
  ([&] { result = (result * kPrime) ^ std::hash<Ts>{}(values); }(), ...);  // NOLINT
  return result;
}

template <typename Iter>
std::uint64_t HashRange(Iter begin, Iter end) {
  using ValueType = typename std::iterator_traits<Iter>::value_type;
  constexpr std::uint64_t kPrime{0x100000001b3};
  std::uint64_t result{0xcbf29ce484222325};
  for (auto iter = begin; iter != end; ++iter) {
    result = (result * kPrime) ^ std::hash<ValueType>{}(*iter);
  }
  return result;
}

}  // namespace xoos::util::hash
