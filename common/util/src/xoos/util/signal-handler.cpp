#include "signal-handler.h"

#include <csignal>
#include <mutex>
#include <ranges>
#include <stdexcept>

namespace xoos::util {

static std::multimap<int, HandlerFunc> handlers;
static std::mutex handlers_mutex;

SignalHandlerHandle AddSignalHandler(int signum, const HandlerFunc& func) {
  const std::lock_guard<std::mutex> lock(handlers_mutex);

  // If this is the first handler we are adding for this signal, set the OS
  // signal handler to call this handler and any future added handlers.
  if (handlers.count(signum) == 0) {
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
