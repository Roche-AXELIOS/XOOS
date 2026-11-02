#include "kstring.h"

#include <string>

#include <htslib/kstring.h>

namespace xoos::io {

Kstring::~Kstring() {
  ks_free(&_ks);
}

Kstring::Kstring(Kstring&& other) noexcept : _ks{other._ks} {
  other._ks = KS_INITIALIZE;
}

Kstring& Kstring::operator=(Kstring&& other) & noexcept {
  if (this != &other) {
    ks_free(&_ks);
    _ks = other._ks;
    other._ks = KS_INITIALIZE;
  }
  return *this;
}

std::string Kstring::Str() {
  return {ks_c_str(&_ks), ks_len(&_ks)};
}

kstring_t* Kstring::Get() {
  return &_ks;
}

void Kstring::Clear() {
  ks_clear(&_ks);
}

}  // namespace xoos::io
