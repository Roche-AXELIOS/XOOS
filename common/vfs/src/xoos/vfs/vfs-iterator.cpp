#include "xoos/vfs/vfs-iterator.h"

#include <cassert>

namespace xoos::vfs {

VirtualFileStream::VirtualFileStream(vfs::VirtualFileHandlePtr handle, ptrdiff_t buffer_size)
    : _handle(std::move(handle)),
      _buffer(AllocateBuffer(buffer_size)),
      _buffer_size{buffer_size},
      _buffer_content_size{0},
      _buffer_offset{0} {
  Advance();
}

void VirtualFileStream::Advance() {
  _buffer_offset++;
  if (_buffer_offset >= _buffer_content_size) {
    _buffer_offset = 0;
    _buffer_content_size = _handle->Read(_buffer.get(), _buffer_size);
    if (_buffer_content_size < 0) {
      throw std::runtime_error(_handle->GetError());
    }
  }
}

const char& VirtualFileStream::GetCurrent() const {
  assert(_buffer_offset < _buffer_size);
  assert(_buffer_offset < _buffer_content_size);
  return _buffer.get()[_buffer_offset];
}

bool VirtualFileStream::Eof() const {
  return _buffer_content_size <= 0;
}

std::unique_ptr<char[]> VirtualFileStream::AllocateBuffer(ptrdiff_t buffer_size) {
  assert(buffer_size > 0);
  return std::unique_ptr<char[]>(new char[buffer_size]);
}

VirtualFileIterator::VirtualFileIterator() : _stream{nullptr} {
}

VirtualFileIterator::VirtualFileIterator(VirtualFileStreamPtr stream) : _stream{std::move(stream)} {
}

VirtualFileIterator& VirtualFileIterator::operator++() {
  assert(_stream != nullptr);
  assert(!_stream->Eof());
  _stream->Advance();
  return *this;
}

const VirtualFileIterator VirtualFileIterator::operator++(int) {  // NOLINT(readability-const-return-type)
  VirtualFileIterator result = *this;
  this->operator++();
  return result;
}

bool VirtualFileIterator::operator!=(const VirtualFileIterator& rhs) const {
  // Handle special case of reaching the end of the file and comparing against the end sentinel iterator
  if (rhs.IsEndSentinel() && _stream->Eof()) {
    return false;
  }
  return rhs._stream != _stream;
}

VirtualFileIterator::reference VirtualFileIterator::operator*() const {
  return _stream->GetCurrent();
}

bool VirtualFileIterator::IsEndSentinel() const {
  return _stream == nullptr;
}

VirtualFileIterator begin(const VirtualFileStreamPtr& stream) {
  return VirtualFileIterator{stream};
}

VirtualFileIterator end(const VirtualFileStreamPtr& stream) {
  return {};
}

VirtualFileIterator begin(const VirtualFileHandlePtr& handle) {
  return VirtualFileIterator{std::make_shared<VirtualFileStream>(handle)};
}

VirtualFileIterator end(const VirtualFileHandlePtr& handle) {
  return {};
}

}  // namespace xoos::vfs
