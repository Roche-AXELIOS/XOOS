#include "xoos/util/tmp-file.h"

#include <filesystem>
#include <fstream>
#include <ios>
#include <stdexcept>
#include <system_error>
#include <utility>

#include <fmt/format.h>

namespace xoos {

namespace fs = std::filesystem;

TmpFile::TmpFile(fs::path path) : _path{std::move(path)} {
  if (!fs::exists(_path)) {
    CreateTmpFile(_path);
  }
}

TmpFile::~TmpFile() {
  Cleanup();
}

const fs::path& TmpFile::Path() const {
  return _path;
}

void TmpFile::CreateTmpFile(const fs::path& path) {
  if (path.has_parent_path()) {
    fs::create_directories(path.parent_path());
  }
  const std::ofstream ofs(path, std::ios::binary);
  if (!ofs) {
    throw std::runtime_error{fmt::format("Could not create temporary file: '{}'", path.string())};
  }
}

void TmpFile::Cleanup() const {
  std::error_code ec;
  fs::remove(_path, ec);
}

}  // namespace xoos
