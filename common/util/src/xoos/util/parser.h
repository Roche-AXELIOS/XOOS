#pragma once
#include <optional>
#include <string>

#include <xoos/types/int.h>

namespace xoos::util {

std::string CreateThreadCountOutOfRangeErrorMessage(const std::string_view& text);

std::string CreateInvalidThreadCountErrorMessage(const std::string_view& text);

std::string CreateNonIntegerThreadCountErrorMessage(
    const std::string_view& text, const std::optional<std::string_view>& detailed_error = std::nullopt);

u64 ParseU64(const std::string& value);

}  // namespace xoos::util
