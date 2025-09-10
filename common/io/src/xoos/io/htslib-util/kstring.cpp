#include "kstring.h"

namespace xoos::io {

Kstring::~Kstring() {
  ks_free(&_ks);
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
