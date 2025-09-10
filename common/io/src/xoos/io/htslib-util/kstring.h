#pragma once

#include <string>

#include <htslib/kstring.h>

namespace xoos::io {

class Kstring {
 public:
  ~Kstring();

  std::string Str();

  kstring_t* Get();

  void Clear();

 private:
  kstring_t _ks = KS_INITIALIZE;
};

}  // namespace xoos::io
