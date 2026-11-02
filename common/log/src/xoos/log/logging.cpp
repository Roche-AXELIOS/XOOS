#include "xoos/log/logging.h"

#include <chrono>
#include <exception>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <spdlog/logger.h>
#include <spdlog/sinks/sink.h>
#include <spdlog/spdlog.h>

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace xoos {

namespace log {

namespace {

template <typename K, typename V>
std::optional<V> Find(const std::unordered_map<K, V>& elements, const K& k) {
  auto it = elements.find(k);
  if (it == elements.end()) {
    return std::nullopt;
  }
  return it->second;
}

}  // namespace

std::optional<LogLevel> ParseLogLevel(const std::string& level) {
  static const auto kLevels = std::unordered_map<std::string, LogLevel>{
      {"debug", LogLevel::kDebug},
      {"info", LogLevel::kInfo},
      {"warn", LogLevel::kWarn},
      {"error", LogLevel::kError},
  };
  return Find(kLevels, level);
}

}  // namespace log

std::unique_ptr<spdlog::logger> Logging::logger;

const std::unordered_map<Logging::LogLevel, Logging::SpdLogLevel> Logging::kLevels = {
    {LogLevel::kDebug, SpdLogLevel::debug},
    {LogLevel::kInfo, SpdLogLevel::info},
    {LogLevel::kWarn, SpdLogLevel::warn},
    {LogLevel::kError, SpdLogLevel::err},
};

Logging::SpdLogLevel Logging::ConvertLogLevel(LogLevel level) {
  return kLevels.at(level);
}

void Logging::Initialize(const std::optional<std::string>& out_file) {
  if (logger != nullptr) {
    return;
  }
  auto sinks = std::vector<std::shared_ptr<spdlog::sinks::sink>>();

  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  sinks.emplace_back(console_sink);

  if (out_file.has_value()) {
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(out_file.value(), true);
    sinks.emplace_back(file_sink);
  }

  logger = std::make_unique<spdlog::logger>(spdlog::logger("multi_sink", sinks.begin(), sinks.end()));
  // RFC 2822 compliant (ISO 8601 has a colon in -07:00), with log level and thread id.
  // Example: [2026-03-20T21:55:33.034-0700] [I] [thread 12345] This is a log message.
  logger->set_pattern("[%Y-%m-%dT%H:%M:%S.%e%z] [%^%L%$] [thread %t] %v");

  static constexpr auto kFlushInterval = std::chrono::seconds(5);
  spdlog::flush_every(kFlushInterval);
}

void Logging::Flush() {
  logger->flush();
}

void Logging::SetLevel(LogLevel level) {
  logger->set_level(ConvertLogLevel(level));
}

/**
 * @brief Get the logger's current log level.
 * @return Log level as an enum value.
 */
log::LogLevel Logging::GetLevel() {
  using enum log::LogLevel;
  switch (logger->level()) {
    case SpdLogLevel::debug:
      return kDebug;
    case SpdLogLevel::info:
      return kInfo;
    case SpdLogLevel::warn:
      return kWarn;
    case SpdLogLevel::err:
      return kError;
    default:
      throw std::runtime_error("Unsupported log level: " + std::to_string(logger->level()));
  }
}

void Logging::Error(const std::string& fmt) noexcept {
  if (logger && logger->level() <= SpdLogLevel::err) {
    logger->error(fmt);
  }
}

bool Logging::ErrorNoExcept(const std::string& fmt) noexcept {
  try {
    Error(fmt);
    return true;
  } catch (...) {
    return false;
  }
}

void Logging::Error(const std::exception& e) noexcept {
  Error(e.what());
}

}  // namespace xoos
