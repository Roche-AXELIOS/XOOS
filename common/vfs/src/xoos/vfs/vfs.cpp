#include "xoos/vfs/vfs.h"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <vector>

#include "xoos/vfs/archive.h"
#include "xoos/vfs/std.h"

namespace xoos::vfs {

ptrdiff_t VirtualFileHandle::Read(std::vector<char>& buffer) {
  // Use new char[] (not make_unique or vector) to avoid zero-initializing the entire capacity;
  // only the bytes actually read are used.
  // NOLINTNEXTLINE(*-avoid-c-arrays, cppcoreguidelines-owning-memory)
  const std::unique_ptr<char[]> tmp_buffer{new char[buffer.capacity()]};
  auto* tmp_buffer_begin = tmp_buffer.get();
  const ptrdiff_t amount_read = Read(tmp_buffer_begin, static_cast<ptrdiff_t>(buffer.capacity()));
  if (amount_read <= 0) {
    return amount_read;
  }
  buffer.assign(tmp_buffer_begin, tmp_buffer_begin + amount_read);
  return amount_read;
}

VirtualFilesystemPtr Open(const fs::path& path) {
  if (is_directory(path)) {
    return std::make_shared<StdVirtualFilesystem>(path);
  }
  return std::make_shared<ArchiveVirtualFilesystem>(path);
}

std::vector<fs::path> VirtualFilesystem::ListRecursive(size_t max_results) const {
  return ListRecursiveImpl(max_results);
}

ptrdiff_t Read(const VirtualFileHandlePtr& handle, char* buffer, ptrdiff_t buffer_size) {
  return handle->Read(buffer, buffer_size);
}

}  // namespace xoos::vfs
