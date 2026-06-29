#pragma once

#include <xoos/types/int.h>

#include <cstring>
#include <fstream>
#include <utility>

namespace xoos::demux {

/**
 * The MemoryFileStream is a simple wrapper around ifstream for efficient
 * random access reading by maintaing an in-memory cache buffer. The cache size
 * can be specified throught the template parameter (N - default 1MB).
 * This class is not thread safe.
 */
template <size_t N = 1024 * 1024>
class MemoryFileStream {
 public:
  explicit MemoryFileStream(std::ifstream& ifs) : _file(ifs) {}

  template <class T>
  void Read(T* ret, u32 size = 1) {
    auto num_bytes = sizeof(T) * size;
    if (_buf_end - _buf_offset >= num_bytes) {
      std::memcpy(ret, _buffer + _buf_offset, num_bytes);
      _buf_offset += num_bytes;
    } else {
      auto num_bytes_1 = _buf_end - _buf_offset;
      if (num_bytes_1 > 0) {
        std::memcpy(ret, _buffer + _buf_offset, num_bytes_1);
      }

      // read next block from file
      if (_file_offset >= _file_size) {
        throw std::runtime_error("not enough bytes to read");
      }

      _file.seekg(_file_offset, std::ios::beg);
      _file.read(reinterpret_cast<char*>(_buffer), std::min(N, static_cast<size_t>(_file_size - _file_offset)));
      auto num_bytes_read = _file.gcount();
      auto num_bytes_2 = num_bytes - num_bytes_1;
      if (std::cmp_less_equal(num_bytes_2, num_bytes_read)) {
        // FIX: The original `ret + num_bytes_1` used pointer arithmetic on T*, which advances
        // by num_bytes_1 * sizeof(T) bytes — not num_bytes_1 bytes.  When sizeof(T) > 1
        // (e.g. Read<u32>), this writes past the destination object, corrupting memory.
        // Cast to char* so the offset is in bytes, matching num_bytes_1's unit.
        std::memcpy(reinterpret_cast<char*>(ret) + num_bytes_1, _buffer, num_bytes_2);
        _file_offset += num_bytes_read;
        _buf_offset = num_bytes_2;
        _buf_end = num_bytes_read;
      } else {
        throw std::runtime_error("not enough bytes to read");
      }
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
  }

  void Skip(size_t size) {
    _buf_offset += size;
    if (_buf_offset >= _buf_end) {
      _file_offset += _buf_offset - _buf_end;
      _buf_offset = 0;
      _buf_end = 0;
    }
  }

  void SetFileSize(std::streamsize file_size) { _file_size = file_size; }

 private:
  std::ifstream& _file;
  std::streamoff _file_offset{0};
  std::streamsize _file_size{0};
  size_t _buf_offset{0};
  size_t _buf_end{0};
  u8 _buffer[N] = {};
};
}  // namespace xoos::demux
