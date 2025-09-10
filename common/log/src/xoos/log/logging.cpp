#include "xoos/log/logging.h"

#include <chrono>
#include <optional>

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace xoos {

namespace log {

template <typename K, typename V>
static std::optional<V> Find(const std::unordered_map<K, V>& elements, const K& k) {
  auto it = elements.find(k);
  if (it == elements.end()) {
    return std::nullopt;
  }
  return it->second;
}

std::optional<LogLevel> ParseLogLevel(const std::string& level) {
  auto levels = std::unordered_map<std::string, LogLevel>{
      {"debug", LogLevel::kDebug},
      {"info", LogLevel::kInfo},
      {"warn", LogLevel::kWarn},
      {"error", LogLevel::kError},
  };
  return Find(levels, level);
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
  logger->set_pattern("[%H:%M:%S.%e %z] [%^%L%$] [thread %t] %v");

  spdlog::flush_every(std::chrono::seconds(5));
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
