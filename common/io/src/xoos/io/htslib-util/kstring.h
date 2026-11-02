#pragma once

#include <string>

#include <htslib/kstring.h>

namespace xoos::io {

class Kstring {
 public:
  // Owns a heap-allocated kstring_t buffer; copying would double-free on destruction.
  // Movable: transfers buffer ownership and leaves the source in an empty, destructible state.
  Kstring() = default;
  Kstring(const Kstring&) = delete;
  Kstring(Kstring&& other) noexcept;
  ~Kstring();
  Kstring& operator=(const Kstring&) = delete;
  Kstring& operator=(Kstring&& other) & noexcept;

  std::string Str();

  kstring_t* Get();

  void Clear();

 private:
  kstring_t _ks = KS_INITIALIZE;
};

}  // namespace xoos::io
