#include "xoos/vfs/mem.h"

#include <algorithm>
#include <ranges>  // NOLINT
#include <utility>

namespace xoos::vfs {

BufferVirtualFileHandle::BufferVirtualFileHandle(const std::vector<char>& buffer) : _buffer{buffer} {
}

BufferVirtualFileHandle::BufferVirtualFileHandle(const std::string& buffer) {
  _buffer.reserve(buffer.size());
  _buffer.assign(std::cbegin(buffer), std::cend(buffer));
}

std::string BufferVirtualFileHandle::GetError() const {
  return "";
}

ptrdiff_t BufferVirtualFileHandle::Read(char* buffer, ptrdiff_t buffer_size) {
  if (std::cmp_equal(_offset, _buffer.size())) {
    return 0;
  }
  auto amount_read = std::min(buffer_size, static_cast<ptrdiff_t>(_buffer.size()) - _offset);
  auto begin = std::begin(_buffer) + _offset;
  std::copy(begin, begin + amount_read, buffer);
  _offset += amount_read;
  return amount_read;
}

MemVirtualFilesystem::MemVirtualFilesystem(MemVirtualFileEntries entries) : _entries{std::move(entries)} {
}

std::vector<fs::path> MemVirtualFilesystem::ListRecursiveImpl(size_t max_results) const {
  auto results = std::vector<fs::path>();
  std::ranges::copy(_entries | std::views::keys | std::views::take(max_results), std::back_inserter(results));
  return results;
}

VirtualFileHandlePtr MemVirtualFilesystem::Open(const fs::path& path) const {
  auto it = _entries.find(path);
  if (it == _entries.end()) {
    return nullptr;
  }
  return it->second(path);
}

std::string MemVirtualFilesystem::GetName() const {
  return "MemFs";
}

}  // namespace xoos::vfs
