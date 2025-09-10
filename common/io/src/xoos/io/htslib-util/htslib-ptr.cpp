#include "htslib-ptr.h"

namespace xoos::io {

void SamFileDeleter::operator()(samFile* h) {
  if (h != nullptr) {
    hts_close(h);
  }
}

void SamHdrDeleter::operator()(sam_hdr_t* h) {
  if (h != nullptr) {
    sam_hdr_destroy(h);
  }
}

void BamHdrDeleter::operator()(bam_hdr_t* h) {
  if (h != nullptr) {
    bam_hdr_destroy(h);
  }
}

void Bam1Deleter::operator()(bam1_t* b) {
  if (b != nullptr) {
    bam_destroy1(b);
  }
}

void FastaIdxDeleter::operator()(faidx_t* f) {
  if (f != nullptr) {
    fai_destroy(f);
  }
}

void HtsIdxDeleter::operator()(hts_idx_t* i) {
  if (i != nullptr) {
    hts_idx_destroy(i);
  }
}

void HtsFileDeleter::operator()(htsFile* f) {
  if (f != nullptr) {
    hts_close(f);
  }
}

HtsRegListDeleter::HtsRegListDeleter(const int num_i) : num_i(num_i) {
}

void HtsRegListDeleter::operator()(gsl::owner<hts_reglist_t*> reglist) const {
  if (reglist != nullptr) {
    hts_reglist_free(reglist, num_i);
  }
}

void HtsItrDeleter::operator()(hts_itr_t* itr) {
  if (itr != nullptr) {
    hts_itr_destroy(itr);
  }
}

void HtsItrMultiDeleter::operator()(hts_itr_multi_t* itr) {
  if (itr != nullptr) {
    hts_itr_multi_destroy(itr);
  }
}

}  // namespace xoos::io
