#pragma once

#include <concepts>
#include <optional>

namespace xoos {

template <typename T>
concept Addable = requires(T a, T b) {
  { a + b } -> std::convertible_to<T>;
};  // NOLINT readability/braces

template <typename T>
concept Incrementable = requires(T a) {
  { ++a } -> std::same_as<T&>;
};  // NOLINT readability/braces

template <typename T>
concept Decrementable = requires(T a) {
  { --a } -> std::same_as<T&>;
};  // NOLINT readability/braces

template <typename T>
  requires Addable<T>
inline std::optional<T> operator+(const std::optional<T>& lhs, const std::optional<T>& rhs) {
  if (!lhs.has_value() && !rhs.has_value()) {
    return std::nullopt;
  }
  if (!lhs.has_value()) {
    return rhs;
  }
  if (!rhs.has_value()) {
    return lhs;
  }
  return *lhs + *rhs;
}

template <typename T>
  requires Addable<T>
inline std::optional<T>& operator+=(std::optional<T>& lhs, const std::optional<T>& rhs) {
  if (lhs.has_value() && rhs.has_value()) {
    lhs = *lhs + *rhs;
  } else if (!lhs.has_value() && rhs.has_value()) {
    lhs = rhs;
  }
  return lhs;
}

template <typename T>
  requires Incrementable<T>
inline std::optional<T>& operator++(std::optional<T>& value) {
  if (value.has_value()) {
    ++(*value);
  }
  return value;
}

template <typename T>
  requires Decrementable<T>
inline std::optional<T>& operator--(std::optional<T>& value) {
  if (value.has_value()) {
    --(*value);
  }
  return value;
}

}  // namespace xoos
