#pragma once

#include <fmt/base.h>

#include "xoos/error/error.h"
#include "xoos/log/logging.h"

namespace xoos::svc {

/**
 * @brief External flag to control whether to treat warn messages as errors.
 * If set to true, any warn messages will throw a runtime error.
 */
extern bool warn_as_error;

/**
 * @brief Same as `Logging::Warn`, but throw a runtime error if warn messages are treated as errors.
 * @param fmt Format string for the warn message.
 * @param args Arguments to format the warn message.
 */
template <typename... Args>
void WarnAsErrorIfSet(fmt::format_string<Args...> fmt, Args&&... args) {
  if (warn_as_error && Logging::GetLevel() <= log::LogLevel::kWarn) {
    throw error::Error(fmt, std::forward<Args>(args)...);
  }
  Logging::Warn(fmt, std::forward<Args>(args)...);
}

}  // namespace xoos::svc
