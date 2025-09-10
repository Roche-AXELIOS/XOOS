#include "xoos/util/tmp-dir.h"

#include <unistd.h>

#include <csignal>
#include <filesystem>

#include <fmt/format.h>

namespace xoos {

namespace fs = std::filesystem;

TmpDir::TmpDir(std::source_location location, bool read_only)
    : _value{CreateTmpDir(fs::path{location.file_name()}.filename())},
      _read_only(read_only),
      _cleanup_handlers{AddCleanupHandler()} {
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

fs::path TmpDir::CreateTmpDir(const std::string& name) {
  auto tmp_dir_template = (fs::temp_directory_path() / (name + ".XXXXXX")).generic_string();
  errno = 0;
  auto* tmp_dir = mkdtemp(tmp_dir_template.data());
  if (errno != 0) {
    throw std::runtime_error{fmt::format("Could not create temporary directory: '{}'", strerror(errno))};
  }
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
  return util::AddSignalHandler({SIGSEGV, SIGABRT, SIGINT, SIGFPE}, [this](int) { Cleanup(); });  // NOLINT
}

}  // namespace xoos
