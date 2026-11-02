#pragma once

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <xoos/types/int.h>

namespace xoos::demux {

/**
 * Buffered sequential reader over an ifstream.
 *
 * Maintains a fixed-size in-memory buffer (default 1 MB) and refills it
 * from the underlying ifstream as needed.  Not thread-safe.
 *
 * Two construction modes:
 *  - Shared fd: references an external ifstream.  A seekg() is issued on
 *               every buffer refill because other readers may have moved
 *               the file position.
 *  - Owned fd:  opens a private ifstream.  No seekg() on refill because
 *               the file position is guaranteed correct — nothing else
 *               touches this fd.  This keeps the kernel's readahead
 *               heuristic in sequential mode, which matters on network
 *               filesystems like Lustre.
 */
template <size_t N = 1024 * 1024>
class MemoryFileStream {
 public:
  /// Shared-fd constructor — caller owns the ifstream lifetime.
  explicit MemoryFileStream(std::ifstream& ifs) : _file(ifs) {
  }

  /// Owned-fd constructor — opens a private ifstream on @p path.
  explicit MemoryFileStream(const std::string& path)
      : _owned_file(std::make_unique<std::ifstream>()), _file(*_owned_file) {
    _file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      _file.open(path, std::ios::binary);
    } catch (const std::ios_base::failure& e) {
      throw std::runtime_error("MemoryFileStream: cannot open " + path + ": " + e.what());
    }
  }

  MemoryFileStream(const MemoryFileStream&) = delete;
  MemoryFileStream(MemoryFileStream&&) noexcept = default;
  MemoryFileStream& operator=(const MemoryFileStream&) = delete;
  MemoryFileStream& operator=(MemoryFileStream&&) = delete;

  template <class T>
  void Read(T* ret, u32 size = 1) {
    auto num_bytes = sizeof(T) * size;
    if (_buf_end - _buf_offset >= num_bytes) {
      std::memcpy(ret, _buffer.data() + _buf_offset, num_bytes);
      _buf_offset += num_bytes;
      return;
    }

    // Drain remaining buffered bytes.
    auto num_bytes_1 = _buf_end - _buf_offset;
    if (num_bytes_1 > 0) {
      std::memcpy(ret, _buffer.data() + _buf_offset, num_bytes_1);
    }

    if (_file_offset >= _file_size) {
      throw std::runtime_error("not enough bytes to read");
    }

    RefillBuffer();

    auto num_bytes_2 = num_bytes - num_bytes_1;
    const auto num_bytes_read = _file.gcount();
    if (std::cmp_less_equal(num_bytes_2, num_bytes_read)) {
      // Byte-offset into ret: char* can alias any type per [basic.lval]/11.
      auto* dest = reinterpret_cast<char*>(ret) + num_bytes_1;  // NOSONAR: S5356 — char aliasing is well-defined
      std::memcpy(dest, _buffer.data(), num_bytes_2);
      _buf_offset = num_bytes_2;
      _buf_end = num_bytes_read;
    } else {
      throw std::runtime_error("not enough bytes to read");
    }
  }

  template <class T>
  T Read() {
    T value;
    Read(&value, 1);
    return value;
  }

  void Reset(std::streamoff new_file_offset) {
    _file_offset = new_file_offset;
    _buf_offset = 0;
    _buf_end = 0;
    // Position the fd so the next refill reads from the right place.
    _file.seekg(_file_offset, std::ios::beg);
  }

  void Skip(size_t size) {
    _buf_offset += size;
    if (_buf_offset >= _buf_end) {
      _file_offset += _buf_offset - _buf_end;
      _buf_offset = 0;
      _buf_end = 0;
    }
  }

  void SetFileSize(const std::streamsize file_size) {
    _file_size = file_size;
  }

 private:
  /// Refill the buffer from the file.  Owned fds skip the seek because
  /// the file position is already correct; shared fds must seek because
  /// another reader may have moved the position.
  void RefillBuffer() {
    if (!_owned_file) {
      _file.seekg(_file_offset, std::ios::beg);
    }
    const auto bytes_remaining = static_cast<size_t>(_file_size - _file_offset);
    _file.read(_buffer.data(), std::min(N, bytes_remaining));
    _file_offset += _file.gcount();
  }

  std::unique_ptr<std::ifstream> _owned_file;
  std::ifstream& _file;
  std::streamoff _file_offset{0};
  std::streamsize _file_size{0};
  size_t _buf_offset{0};
  size_t _buf_end{0};
  std::array<char, N> _buffer{};
};
}  // namespace xoos::demux
