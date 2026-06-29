#pragma once

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "vfs.h"

namespace xoos::vfs {

class BufferVirtualFileHandle : public vfs::VirtualFileHandle {
 public:
  explicit BufferVirtualFileHandle(const std::vector<char>& buffer);

  explicit BufferVirtualFileHandle(const std::string& buffer);

  std::string GetError() const override;

  ptrdiff_t Read(char* buffer, ptrdiff_t buffer_size) override;

 private:
  std::vector<char> _buffer;
  ptrdiff_t _offset = 0;
};

class MemVirtualFilesystem : public VirtualFilesystem {
 public:
  using MemVirtualFileEntries = std::map<fs::path, std::function<VirtualFileHandlePtr(const fs::path&)>>;

  explicit MemVirtualFilesystem(MemVirtualFileEntries entries);

  ~MemVirtualFilesystem() override = default;

  std::vector<fs::path> ListRecursiveImpl(size_t max_results) const override;

  VirtualFileHandlePtr Open(const fs::path& path) const override;

  std::string GetName() const override;

 private:
  MemVirtualFileEntries _entries;
};

}  // namespace xoos::vfs
