#pragma once

#include <filesystem>
#include <source_location>  // NOLINT
#include <string>

#include "signal-handler.h"

namespace xoos {

class TmpDir {
 public:
  explicit TmpDir(std::source_location location = std::source_location::current(), bool read_only = false);

  ~TmpDir();

  const std::filesystem::path& Path() const;

  template <class T>
  std::filesystem::path operator/(T&& other) const {
    return _value / other;
  }

 private:
  static std::filesystem::path CreateTmpDir(const std::string& name);

  std::vector<util::SignalHandlerHandle> AddCleanupHandler();
  void Cleanup();

  std::filesystem::path _value;
  bool _read_only;
  std::vector<util::SignalHandlerHandle> _cleanup_handlers;
};

}  // namespace xoos
