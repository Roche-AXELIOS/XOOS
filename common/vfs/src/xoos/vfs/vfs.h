#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace xoos::vfs {

namespace fs = std::filesystem;

class VirtualFileHandle {
 public:
  virtual ~VirtualFileHandle() = default;

  virtual std::string GetError() const = 0;

  /**
   * A standard C-like read API. Will read at most buffer_size char into buffer,
   * but may not fill buffer completely. Return value indicates number of char read.
   * @param buffer the buffer to read into
   * @param buffer_size  the size of the buffer
   * @return -1 indicates error, 0 indicates EOF, >0 is number of char read
   */
  virtual ptrdiff_t Read(char* buffer, ptrdiff_t buffer_size) = 0;

  template <size_t Nm>
  ptrdiff_t Read(const std::array<char, Nm>& buffer) {
    return Read(buffer, buffer.size());
  }

  ptrdiff_t Read(std::vector<char>& buffer);
};

using VirtualFileHandlePtr = std::shared_ptr<VirtualFileHandle>;

class VirtualFilesystem {
 public:
  virtual ~VirtualFilesystem() = default;

  virtual VirtualFileHandlePtr Open(const fs::path& path) const = 0;

  std::vector<fs::path> ListRecursive(size_t max_results = 1000) const;

  virtual std::string GetName() const = 0;

 protected:
  /**
   * Recursively list the contents of the filesystem, directories have a trailing slash.
   * @param max_results the maximum number of results to return, this is a simple precaution to avoid any infinite loops
   */
  virtual std::vector<fs::path> ListRecursiveImpl(size_t max_results) const = 0;
};

using VirtualFilesystemPtr = std::shared_ptr<VirtualFilesystem>;

VirtualFilesystemPtr Open(const fs::path& path = fs::current_path());

ptrdiff_t Read(const VirtualFileHandlePtr& handle, char* buffer, ptrdiff_t buffer_size);

}  // namespace xoos::vfs
