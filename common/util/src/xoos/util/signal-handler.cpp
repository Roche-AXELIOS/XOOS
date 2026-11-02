#include "signal-handler.h"

#include <csignal>
#include <initializer_list>
#include <map>
#include <mutex>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace xoos::util {

namespace {

// Process-wide registry of signal handlers, guarded by handlers_mutex.
// Both must be mutable file-scope state: the C std::signal callback cannot
// capture, so the dispatcher relies on this global registry. NOLINT covers
// clang-tidy; NOSONAR covers the SonarQube non-const-global rule.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::multimap<int, HandlerFunc> handlers;  // NOSONAR
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex handlers_mutex;  // NOSONAR

}  // namespace

SignalHandlerHandle AddSignalHandler(int signum, const HandlerFunc& func) {
  const std::lock_guard<std::mutex> lock(handlers_mutex);

  // If this is the first handler we are adding for this signal, set the OS
  // signal handler to call this handler and any future added handlers.
  if (!handlers.contains(signum)) {
    auto result = std::signal(signum, [](int sig) {
      const std::lock_guard<std::mutex> lock(handlers_mutex);
      auto range = handlers.equal_range(sig);

      // Call in reverse registration order, so that newer handlers are called first. Like unwinding a stack.
      for (const auto& [key, value] : std::ranges::subrange(range.first, range.second) | std::views::reverse) {
        value(sig);
      }
    });

    if (result == SIG_ERR) {
      throw std::runtime_error("Failed to set signal handler");
    }
    handlers.emplace(signum, result);
  }

  return handlers.emplace(signum, func);
}

std::vector<SignalHandlerHandle> AddSignalHandler(std::initializer_list<int> signums, const HandlerFunc& func) {
  std::vector<SignalHandlerHandle> handles;
  for (const auto& signum : signums) {
    handles.push_back(AddSignalHandler(signum, func));
  }
  return handles;
}

void RemoveSignalHandler(const SignalHandlerHandle& handle) {
  const std::lock_guard<std::mutex> lock(handlers_mutex);
  handlers.erase(handle);
}

void RemoveSignalHandler(const std::vector<SignalHandlerHandle>& handles) {
  for (const auto& handle : handles) {
    RemoveSignalHandler(handle);
  }
}

}  // namespace xoos::util
