#include "htslib-bcf-ptr.h"

#include <htslib/tbx.h>
#include <htslib/vcf.h>

namespace xoos::io {

void BcfHdrDeleter::operator()(bcf_hdr_t* const h) const {
  if (h != nullptr) {
    bcf_hdr_destroy(h);
  }
}

void Bcf1Deleter::operator()(bcf1_t* const b) const {
  if (b != nullptr) {
    bcf_destroy(b);
  }
}

void TbxDeleter::operator()(tbx_t* const t) const {
  if (t != nullptr) {
    tbx_destroy(t);
  }
}

}  // namespace xoos::io
