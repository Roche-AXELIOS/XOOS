#include "xoos/util/tmp-dir.h"

#include <cerrno>
#include <csignal>
#include <cstring>
#include <filesystem>
#include <source_location>
#include <stdexcept>
#include <string>
#include <vector>

#include <fmt/format.h>

#include "signal-handler.h"

namespace xoos {

namespace fs = std::filesystem;

TmpDir::TmpDir(const std::source_location location, const bool read_only)
    : _value{CreateTmpDir(fs::temp_directory_path(), fs::path{location.file_name()}.filename())},
      _read_only(read_only),
      _cleanup_handlers{AddCleanupHandler()} {
  ApplyReadOnlyPermissions();
}

TmpDir::TmpDir(const fs::path& parent_dir, const std::string& prefix) : TmpDir(parent_dir, prefix, false) {
}

TmpDir::TmpDir(const fs::path& parent_dir, const std::string& prefix, const bool read_only)
    : _value{CreateTmpDir(parent_dir, prefix)}, _read_only(read_only), _cleanup_handlers{AddCleanupHandler()} {
  ApplyReadOnlyPermissions();
}

void TmpDir::ApplyReadOnlyPermissions() const {
  if (_read_only) {
    try {
      fs::permissions(
          _value, fs::perms::owner_write | fs::perms::group_write | fs::perms::others_write, fs::perm_options::remove);
    } catch (const fs::filesystem_error& e) {
      throw std::runtime_error{fmt::format("Unable to make temporary directory read-only: '{}'", e.what())};
    }
  }
}

TmpDir::~TmpDir() {
  Cleanup();
  util::RemoveSignalHandler(_cleanup_handlers);
}

const fs::path& TmpDir::Path() const {
  return _value;
}

fs::path TmpDir::CreateTmpDir(const fs::path& parent_dir, const std::string& prefix) {
  auto tmp_dir_template = (parent_dir / (prefix + ".XXXXXX")).string();
  // NOLINTNEXTLINE(misc-include-cleaner) - mkdtemp is the POSIX form declared via <cstdlib>.
  const auto* const tmp_dir = mkdtemp(tmp_dir_template.data());
  if (tmp_dir == nullptr) {
    // NOLINTNEXTLINE(concurrency-mt-unsafe) - strerror is called only on the fatal mkdtemp failure path.
    throw std::runtime_error{fmt::format("Could not create temporary directory: '{}'", strerror(errno))};
  }
  // mkdtemp creates the directory with 0700 (owner-only) permissions, which is safe
  // even under a publicly writable parent like /tmp. Enforce explicitly for defense-in-depth.
  fs::permissions(fs::path{tmp_dir}, fs::perms::owner_all, fs::perm_options::replace);
  return fs::path{tmp_dir};
}

void TmpDir::Cleanup() {
  if (_read_only) {
    try {
      fs::permissions(_value, fs::perms::owner_write, fs::perm_options::add);
    } catch (const fs::filesystem_error& e) {
      throw std::runtime_error{fmt::format("Unable to make temporary directory writable: '{}'", e.what())};
    }
  }
  fs::remove_all(_value);
}

std::vector<util::SignalHandlerHandle> TmpDir::AddCleanupHandler() {
  return util::AddSignalHandler(  // NOLINT
      {SIGSEGV, SIGABRT, SIGINT, SIGFPE},
      [this](const int /*signum*/) {  // NOSONAR(M23_058, S813)
        Cleanup();
      });
}

}  // namespace xoos
