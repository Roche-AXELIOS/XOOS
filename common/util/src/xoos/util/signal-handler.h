#pragma once

#include <functional>
#include <map>
#include <vector>

namespace xoos::util {

using HandlerFunc = std::function<void(int)>;
using SignalHandlerHandle = std::multimap<int, HandlerFunc>::const_iterator;

SignalHandlerHandle AddSignalHandler(int signum, const HandlerFunc& func);

std::vector<SignalHandlerHandle> AddSignalHandler(std::initializer_list<int> signums, const HandlerFunc& func);

void RemoveSignalHandler(const SignalHandlerHandle& handle);

void RemoveSignalHandler(const std::vector<SignalHandlerHandle>& handles);

}  // namespace xoos::util
