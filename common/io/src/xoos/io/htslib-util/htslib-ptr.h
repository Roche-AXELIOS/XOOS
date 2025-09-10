#pragma once

#include <gsl/gsl>
#include <memory>

#include <htslib/faidx.h>
#include <htslib/sam.h>

namespace xoos::io {

struct SamFileDeleter {
  void operator()(samFile* h);
};

using SamFilePtr = std::unique_ptr<samFile, SamFileDeleter>;

struct SamHdrDeleter {
  void operator()(sam_hdr_t* h);
};

using SamHdrPtr = std::unique_ptr<sam_hdr_t, SamHdrDeleter>;

struct BamHdrDeleter {
  void operator()(bam_hdr_t* h);
};

using BamHdrPtr = std::unique_ptr<bam_hdr_t, BamHdrDeleter>;

struct Bam1Deleter {
  void operator()(bam1_t* b);
};

using Bam1Ptr = std::unique_ptr<bam1_t, Bam1Deleter>;

struct FastaIdxDeleter {
  void operator()(faidx_t* f);
};

using FaIdxPtr = std::unique_ptr<faidx_t, FastaIdxDeleter>;

struct HtsIdxDeleter {
  void operator()(hts_idx_t* i);
};

using HtsIdxPtr = std::unique_ptr<hts_idx_t, HtsIdxDeleter>;

struct HtsItrDeleter {
  void operator()(hts_itr_t* itr);
};

using HtsItrPtr = std::unique_ptr<hts_itr_t, HtsItrDeleter>;

struct HtsItrMultiDeleter {
  void operator()(hts_itr_multi_t* itr);
};

using HtsItrMultiPtr = std::unique_ptr<hts_itr_multi_t, HtsItrMultiDeleter>;

struct HtsFileDeleter {
  void operator()(htsFile* f);
};

using HtsFilePtr = std::unique_ptr<htsFile, HtsFileDeleter>;

struct HtsRegListDeleter {
  // number of elements to delete
  // we track so we can free the reglist correctly
  int num_i{};

  explicit HtsRegListDeleter(int num_i);

  void operator()(gsl::owner<hts_reglist_t*> reglist) const;
};

using HtsRegListPtr = std::unique_ptr<hts_reglist_t, HtsRegListDeleter>;

}  // namespace xoos::io
