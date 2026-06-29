#pragma once

#include <memory>
#include <string>
#include <vector>

#include "vfs.h"

struct archive;

namespace xoos::vfs {

class ArchiveVirtualFileHandle : public VirtualFileHandle {
 public:
  explicit ArchiveVirtualFileHandle(const std::shared_ptr<struct archive>& archive);

  ~ArchiveVirtualFileHandle() override = default;

  std::string GetError() const override;

  ptrdiff_t Read(char* buffer, ptrdiff_t buffer_size) override;

 private:
  std::shared_ptr<struct archive> _archive;
};

class ArchiveVirtualFilesystem : public VirtualFilesystem {
 public:
  explicit ArchiveVirtualFilesystem(fs::path archive, size_t block_size = 10240);

  ~ArchiveVirtualFilesystem() override = default;

  std::vector<fs::path> ListRecursiveImpl(size_t max_results) const override;

  VirtualFileHandlePtr Open(const fs::path& path) const override;

  std::string GetName() const override;

 private:
  using ArchivePtr = std::shared_ptr<struct archive>;

  ArchivePtr OpenArchive() const;

  static std::string GetError(const ArchivePtr& archive);

 private:
  fs::path _archive;
  size_t _block_size;
};

}  // namespace xoos::vfs
