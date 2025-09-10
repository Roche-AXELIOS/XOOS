#pragma once

#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>

#include <spdlog/spdlog.h>

#define XOOS_LOG_ONCE_FLAG_VARNAME_CONCAT(prefix, line) prefix##line
#define XOOS_LOG_ONCE_FLAG_VARNAME(prefix, line) XOOS_LOG_ONCE_FLAG_VARNAME_CONCAT(prefix, line)

// Same as `Logging::Debug`, but only logs the message once. Repeated calls from the same call site will not log the
// message.
#define XOOS_LOG_DEBUG_ONCE(...)                                                        \
  static std::once_flag XOOS_LOG_ONCE_FLAG_VARNAME(__xoos_log_debug_logged_, __LINE__); \
  std::call_once(XOOS_LOG_ONCE_FLAG_VARNAME(__xoos_log_debug_logged_, __LINE__), [&] { Logging::Debug(__VA_ARGS__); });

// Same as `Logging::Info`, but only logs the message once. Repeated calls from the same call site will not log the
// message.
#define XOOS_LOG_INFO_ONCE(...)                                                        \
  static std::once_flag XOOS_LOG_ONCE_FLAG_VARNAME(__xoos_log_info_logged_, __LINE__); \
  std::call_once(XOOS_LOG_ONCE_FLAG_VARNAME(__xoos_log_info_logged_, __LINE__), [&] { Logging::Info(__VA_ARGS__); });

// Same as `Logging::Warn`, but only logs the message once. Repeated calls from the same call site will not log the
// message.
#define XOOS_LOG_WARN_ONCE(...)                                                        \
  static std::once_flag XOOS_LOG_ONCE_FLAG_VARNAME(__xoos_log_warn_logged_, __LINE__); \
  std::call_once(XOOS_LOG_ONCE_FLAG_VARNAME(__xoos_log_warn_logged_, __LINE__), [&] { Logging::Warn(__VA_ARGS__); });

// Same as `Logging::Error`, but only logs the message once. Repeated calls from the same call site will not log the
// message.
#define XOOS_LOG_ERROR_ONCE(...)                                                        \
  static std::once_flag XOOS_LOG_ONCE_FLAG_VARNAME(__xoos_log_error_logged_, __LINE__); \
  std::call_once(XOOS_LOG_ONCE_FLAG_VARNAME(__xoos_log_error_logged_, __LINE__), [&] { Logging::Error(__VA_ARGS__); });

namespace xoos {

namespace log {

enum class LogLevel {
  kDebug,
  kInfo,
  kWarn,
  kError
};

std::optional<LogLevel> ParseLogLevel(const std::string& level);

}  // namespace log

class Logging {
 private:  // NOLINT - private access is intentional to provide a singleton-like interface below for template arguments.
  static std::unique_ptr<spdlog::logger> logger;

  using LogLevel = log::LogLevel;
  using SpdLogLevel = spdlog::level::level_enum;

  static const std::unordered_map<LogLevel, SpdLogLevel> kLevels;

  static SpdLogLevel ConvertLogLevel(LogLevel level);

 public:
  Logging() = delete;

  static void Initialize(const std::optional<std::string>& out_file = std::nullopt);
  static void Flush();
  static LogLevel GetLevel();
  static void SetLevel(LogLevel level);

  template <typename... Args>
  static void Log(LogLevel level, fmt::format_string<Args...> fmt, Args&&... args) {
    if (logger) {
      logger->log(ConvertLogLevel(level), fmt, std::forward<Args>(args)...);
    }
  }

  template <typename... Args>
  static void Debug(fmt::format_string<Args...> fmt, Args&&... args) {
    // The logger builds an async queue, which will begin to block when it becomes full.
    // Additionally, we want to avoid calling the format string call for strings that are never logged.
    // This will be the same for every log level.
    if (logger && logger->level() <= SpdLogLevel::debug) {
      logger->debug(fmt, std::forward<Args>(args)...);
    }
  }

  template <typename... Args>
  static bool DebugNoExcept(fmt::format_string<Args...> fmt, Args&&... args) noexcept {
    try {
      Debug(fmt, std::forward<Args>(args)...);
      return true;
    } catch (...) {
      return false;
    }
  }

  template <typename... Args>
  static void Info(fmt::format_string<Args...> fmt, Args&&... args) {
    if (logger && logger->level() <= SpdLogLevel::info) {
      logger->info(fmt, std::forward<Args>(args)...);
    }
  }

  template <typename... Args>
  static bool InfoNoExcept(fmt::format_string<Args...> fmt, Args&&... args) noexcept {
    try {
      Info(fmt, std::forward<Args>(args)...);
      return true;
    } catch (...) {
      return false;
    }
  }

  template <typename... Args>
  static void Warn(fmt::format_string<Args...> fmt, Args&&... args) {
    if (logger && logger->level() <= SpdLogLevel::warn) {
      logger->warn(fmt, std::forward<Args>(args)...);
    }
  }

  template <typename... Args>
  static bool WarnNoExcept(fmt::format_string<Args...> fmt, Args&&... args) noexcept {
    try {
      Warn(fmt, std::forward<Args>(args)...);
      return true;
    } catch (...) {
      return false;
    }
  }

  template <typename... Args>
  static void Error(fmt::format_string<Args...> fmt, Args&&... args) {
    if (logger && logger->level() <= SpdLogLevel::err) {
      logger->error(fmt, std::forward<Args>(args)...);
    }
  }

  // Log an an error message but do not throw an exception, useful for error handling in destructors or other
  // situations where exceptions are not allowed.
  template <typename... Args>
  static bool ErrorNoExcept(fmt::format_string<Args...> fmt, Args&&... args) noexcept {
    try {
      Error(fmt, std::forward<Args>(args)...);
      return true;
    } catch (...) {
      return false;
    }
  }

  static void Error(const std::string& fmt) noexcept;

  static bool ErrorNoExcept(const std::string& fmt) noexcept;

  static void Error(const std::exception& e) noexcept;
};

}  // namespace xoos
