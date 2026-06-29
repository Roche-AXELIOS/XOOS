#pragma once

#include <memory>

#include "vfs.h"

namespace xoos::vfs {

class VirtualFileStream {
 public:
  explicit VirtualFileStream(vfs::VirtualFileHandlePtr handle, ptrdiff_t buffer_size = 4096);

  void Advance();

  const char& GetCurrent() const;

  bool Eof() const;

 private:
  static std::unique_ptr<char[]> AllocateBuffer(ptrdiff_t buffer_size);

  vfs::VirtualFileHandlePtr _handle;

  std::unique_ptr<char[]> _buffer;
  ptrdiff_t _buffer_size;
  ptrdiff_t _buffer_content_size;
  ptrdiff_t _buffer_offset;
};

using VirtualFileStreamPtr = std::shared_ptr<VirtualFileStream>;

class VirtualFileIterator {
 public:
  // The following types are required to implement the LegacyInputIterator concept
  // NOLINTBEGIN
  using difference_type = std::ptrdiff_t;
  using value_type = char;
  using pointer = std::add_pointer_t<std::add_const_t<value_type>>;
  using reference = std::add_lvalue_reference_t<std::add_const_t<value_type>>;
  using iterator_category = std::input_iterator_tag;
  // NOLINTEND

  VirtualFileIterator();

  explicit VirtualFileIterator(VirtualFileStreamPtr stream);

  VirtualFileIterator& operator++();

  const VirtualFileIterator operator++(int);

  bool operator!=(const VirtualFileIterator& rhs) const;

  reference operator*() const;

 private:
  bool IsEndSentinel() const;

  VirtualFileStreamPtr _stream;
};

// Enforce that VirtualFileIterator implements the LegacyInputIterator concept
static_assert(std::input_iterator<VirtualFileIterator>);

VirtualFileIterator begin(const VirtualFileStreamPtr& stream);  // NOLINT(readability-identifier-naming)
VirtualFileIterator end(const VirtualFileStreamPtr& stream);    // NOLINT(readability-identifier-naming)

VirtualFileIterator begin(const VirtualFileHandlePtr& handle);  // NOLINT(readability-identifier-naming)
VirtualFileIterator end(const VirtualFileHandlePtr& handle);    // NOLINT(readability-identifier-naming)

}  // namespace xoos::vfs
