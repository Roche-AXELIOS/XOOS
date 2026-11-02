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

  ArchiveVirtualFileHandle(const ArchiveVirtualFileHandle&) = default;
  ArchiveVirtualFileHandle(ArchiveVirtualFileHandle&&) noexcept = default;
  ArchiveVirtualFileHandle& operator=(const ArchiveVirtualFileHandle&) = default;
  ArchiveVirtualFileHandle& operator=(ArchiveVirtualFileHandle&&) noexcept = default;
  ~ArchiveVirtualFileHandle() override = default;

  std::string GetError() const override;

  ptrdiff_t Read(char* buffer, ptrdiff_t buffer_size) override;

 private:
  std::shared_ptr<struct archive> _archive;
};

class ArchiveVirtualFilesystem : public VirtualFilesystem {
 public:
  static constexpr size_t kDefaultBlockSize = 10'240;

  explicit ArchiveVirtualFilesystem(fs::path archive, size_t block_size = kDefaultBlockSize);

  ArchiveVirtualFilesystem(const ArchiveVirtualFilesystem&) = default;
  ArchiveVirtualFilesystem(ArchiveVirtualFilesystem&&) noexcept = default;
  ArchiveVirtualFilesystem& operator=(const ArchiveVirtualFilesystem&) = default;
  ArchiveVirtualFilesystem& operator=(ArchiveVirtualFilesystem&&) noexcept = default;
  ~ArchiveVirtualFilesystem() override = default;

  std::vector<fs::path> ListRecursiveImpl(size_t max_results) const override;

  VirtualFileHandlePtr Open(const fs::path& path) const override;

  std::string GetName() const override;

 private:
  using ArchivePtr = std::shared_ptr<struct archive>;

  ArchivePtr OpenArchive() const;

  static std::string GetError(const ArchivePtr& archive);

  fs::path _archive;
  size_t _block_size;
};

}  // namespace xoos::vfs
