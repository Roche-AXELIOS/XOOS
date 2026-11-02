#include "htslib-ptr.h"

#include <gsl/pointers>
#include <memory>

#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/util/adopt.h>

namespace xoos::io {

void SamFileDeleter::operator()(samFile* const h) const {
  if (h != nullptr) {
    hts_close(h);
  }
}

void SamHdrDeleter::operator()(sam_hdr_t* const h) const {
  if (h != nullptr) {
    sam_hdr_destroy(h);
  }
}

void BamHdrDeleter::operator()(bam_hdr_t* const h) const {
  if (h != nullptr) {
    bam_hdr_destroy(h);
  }
}

void Bam1Deleter::operator()(bam1_t* const b) const {
  if (b != nullptr) {
    bam_destroy1(b);
  }
}

void FastaIdxDeleter::operator()(faidx_t* const f) const {
  if (f != nullptr) {
    fai_destroy(f);
  }
}

void HtsIdxDeleter::operator()(hts_idx_t* const i) const {
  if (i != nullptr) {
    hts_idx_destroy(i);
  }
}

void HtsFileDeleter::operator()(htsFile* const f) const {
  if (f != nullptr) {
    hts_close(f);
  }
}

HtsRegListDeleter::HtsRegListDeleter(const s32 num_i) : _num_i(num_i) {
}

void HtsRegListDeleter::operator()(const gsl::owner<hts_reglist_t*> reglist) const {
  if (reglist != nullptr) {
    hts_reglist_free(reglist, _num_i);
  }
}

void HtsItrDeleter::operator()(hts_itr_t* const itr) const {
  if (itr != nullptr) {
    hts_itr_destroy(itr);
  }
}

void HtsItrMultiDeleter::operator()(hts_itr_multi_t* const itr) const {
  if (itr != nullptr) {
    hts_itr_multi_destroy(itr);
  }
}

std::shared_ptr<bam1_t> MakeBam1SharedPtr(bam1_t* const ptr) {
  return AdoptShared(ptr, Bam1Deleter{});
}

// misc-misplaced-const: bam_plp_t is a pointer typedef; const on the typedef correctly
// const-qualifies the pointer value, but clang-tidy incorrectly treats this as misplaced.
// NOLINTBEGIN(misc-misplaced-const)
void BamPlpDeleter::operator()(const bam_plp_t p) const {
  // NOLINTEND(misc-misplaced-const)
  if (p != nullptr) {
    bam_plp_destroy(p);
  }
}

void BgzfDeleter::operator()(BGZF* const f) const {
  if (f != nullptr) {
    bgzf_close(f);
  }
}

}  // namespace xoos::io
