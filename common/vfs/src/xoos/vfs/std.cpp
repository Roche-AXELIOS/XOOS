#include "xoos/vfs/std.h"

#include <cstring>

namespace xoos::vfs {

StdVirtualFileHandle::StdVirtualFileHandle(std::ifstream stream) : _stream{std::move(stream)} {
}

std::string StdVirtualFileHandle::GetError() const {
  return strerror(errno);
}

ptrdiff_t StdVirtualFileHandle::Read(char* buffer, ptrdiff_t buffer_size) {
  const std::streamsize amount_read = _stream.readsome(buffer, buffer_size);
  if (_stream.fail() || _stream.bad()) {
    return -1;
  }
  return amount_read;
}

StdVirtualFilesystem::StdVirtualFilesystem(fs::path root) : _root(std::move(root)) {
}

std::vector<fs::path> StdVirtualFilesystem::ListRecursiveImpl(size_t max_results) const {
  std::vector<fs::path> results;
  for (const auto& p : fs::recursive_directory_iterator{_root}) {
    if (results.size() == max_results) {
      break;
    }
    results.emplace_back(Normalize(p));
  }
  return results;
}

VirtualFileHandlePtr StdVirtualFilesystem::Open(const fs::path& path) const {
  auto full_path = _root / path;
  if (!fs::exists(full_path)) {
    return nullptr;
  }
  return std::make_shared<StdVirtualFileHandle>(std::ifstream{full_path});
}

std::string StdVirtualFilesystem::GetName() const {
  return _root.string();
}

static fs::path AddTrailingSeparator(const fs::path& path) {
  return fs::path{path.generic_string() + fs::path::preferred_separator};
}

fs::path StdVirtualFilesystem::Normalize(const fs::path& path) const {
  auto relative_path = fs::relative(path, _root);
  return fs::is_directory(path) ? AddTrailingSeparator(relative_path) : relative_path;
}

}  // namespace xoos::vfs
