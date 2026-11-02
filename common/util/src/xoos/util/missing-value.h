#pragma once

#include <optional>
#include <string>
#include <string_view>

#include <fmt/format.h>

namespace xoos::util {

/**
 * @brief XOOS standard tokens for values that are missing or not applicable in module output.
 *
 * Each output format has exactly one token, so modules must not define their own. See "Value Representation" in
 * docs/user-experience-guide/command-line-interface-best-practices.md for the convention these constants encode.
 */

/// Missing value in a TSV column or a free-form BED column. The XOOS default.
inline constexpr std::string_view kNA = "NA";

/// Missing value in a spec-governed field: VCF columns and the BED fields defined by the BED specification.
inline constexpr std::string_view kMissingSpecField = ".";

/**
 * @brief Check whether a field read back from an output file is the missing-value token.
 *
 * The comparison is exact. kNA is the only supported spelling: alternatives such as "nan" or "N/A" are deliberately
 * not recognized, so input that predates this convention fails to parse rather than being silently reinterpreted.
 *
 * @param value field to test
 * @return true if the field is exactly kNA
 */
bool IsNA(std::string_view value);  // NOSONAR - exported helper used by downstream modules/tests

/**
 * @brief Format an optional value for TSV/BED output, substituting kNA when it is absent.
 *
 * @param value optional to format
 * @param format_string compile-time fmt format string applied to the contained value, e.g. "{:.2f}"
 * @return the formatted value, or kNA if the optional is empty
 */
template <typename T>
std::string FormatOrNA(const std::optional<T>& value,  // NOSONAR - exported helper used by downstream modules/tests
                       fmt::format_string<const T&> format_string) {
  if (!value.has_value()) {
    return std::string{kNA};
  }
  return fmt::format(format_string, *value);
}

}  // namespace xoos::util
