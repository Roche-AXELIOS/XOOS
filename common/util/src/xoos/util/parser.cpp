#include "parser.h"

#include <algorithm>
#include <optional>

#include <fmt/format.h>

#include <xoos/types/int.h>

namespace xoos::util {

std::string CreateThreadCountOutOfRangeErrorMessage(const std::string_view& text) {
  return fmt::format("'{}' cannot be converted to option because it is out of range of the maximum thread size.", text);
}

std::string CreateInvalidThreadCountErrorMessage(const std::string_view& text) {
  return fmt::format("'{}' cannot be converted to option because it is an invalid maximum thread size.", text);
}

std::string CreateNonIntegerThreadCountErrorMessage(const std::string_view& text,
                                                    const std::optional<std::string_view>& detailed_error) {
  const auto msg = fmt::format("'{}' cannot be converted to option, must be a positive integer.", text);
  const auto detailed_msg = detailed_error ? fmt::format(" Error: \n{}", *detailed_error) : "";
  return fmt::format("{}{}", msg, detailed_msg);
}

u64 ParseU64(const std::string& value) {
  // Check that all characters are digits
  if (!std::ranges::all_of(value, ::isdigit)) {
    throw std::invalid_argument(CreateNonIntegerThreadCountErrorMessage(value));
  }

  try {
    // Convert the string to u64, which is unsigned and should not be negative
    return std::stoull(value);
  } catch (const std::invalid_argument&) {
    throw std::invalid_argument(CreateInvalidThreadCountErrorMessage(value));
  } catch (const std::out_of_range&) {
    throw std::invalid_argument(CreateThreadCountOutOfRangeErrorMessage(value));
  } catch (const std::exception& e) {
    // Catch any other exceptions and rethrow with a more specific message
    throw std::runtime_error(CreateNonIntegerThreadCountErrorMessage(value, fmt::format("Exception: {}", e.what())));
  }
}
}  // namespace xoos::util
