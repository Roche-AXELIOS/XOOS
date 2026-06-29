#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "vfs.h"

namespace xoos::vfs {

class StdVirtualFileHandle : public VirtualFileHandle {
 public:
  explicit StdVirtualFileHandle(std::ifstream stream);

  ~StdVirtualFileHandle() override = default;

  std::string GetError() const override;

  ptrdiff_t Read(char* buffer, ptrdiff_t buffer_size) override;

 private:
  std::ifstream _stream;
};

class StdVirtualFilesystem : public VirtualFilesystem {
 public:
  explicit StdVirtualFilesystem(fs::path root = std::filesystem::current_path());

  ~StdVirtualFilesystem() override = default;

  std::vector<fs::path> ListRecursiveImpl(size_t max_results) const override;

  VirtualFileHandlePtr Open(const fs::path& path) const override;

  std::string GetName() const override;

 private:
  /**
   * Make path relative to root of the filesystem, add trailing slash if directory.
   */
  fs::path Normalize(const fs::path& path) const;

 private:
  fs::path _root;
};

}  // namespace xoos::vfs
