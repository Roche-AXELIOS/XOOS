#pragma once

#include <filesystem>
#include <source_location>  // NOLINT
#include <string>
#include <utility>
#include <vector>

#include "signal-handler.h"

namespace xoos {

class TmpDir {
 public:
  explicit TmpDir(  // NOSONAR — used by subclasses and other modules (S5536)
      std::source_location location = std::source_location::current(),  // NOSONAR — must be a default arg (S1712)
      bool read_only = false);                                          // NOSONAR(S1712)

  /**
   * @brief Create a temp directory under @p parent_dir with the given @p prefix.
   *
   * Uses POSIX mkdtemp for atomic, collision-free directory creation.
   * The directory name will be `<prefix>.XXXXXX` where X is replaced by mkdtemp.
   */
  TmpDir(const std::filesystem::path& parent_dir, const std::string& prefix);                  // NOSONAR(S5536)
  TmpDir(const std::filesystem::path& parent_dir, const std::string& prefix, bool read_only);  // NOSONAR(S5536)

  ~TmpDir();

  // A TmpDir owns the lifetime of its temporary directory, so it is neither copyable nor movable.
  TmpDir(const TmpDir&) = delete;
  TmpDir(TmpDir&&) = delete;
  TmpDir& operator=(const TmpDir&) = delete;
  TmpDir& operator=(TmpDir&&) = delete;

  const std::filesystem::path& Path() const;

  template <class T>
  std::filesystem::path operator/(T&& other) const {
    return _value / std::forward<T>(other);
  }

 private:
  static std::filesystem::path CreateTmpDir(const std::filesystem::path& parent_dir, const std::string& prefix);

  void ApplyReadOnlyPermissions() const;

  std::vector<util::SignalHandlerHandle> AddCleanupHandler();
  void Cleanup();

  std::filesystem::path _value;
  bool _read_only;
  std::vector<util::SignalHandlerHandle> _cleanup_handlers;
};

}  // namespace xoos
