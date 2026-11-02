#pragma once

#include <gsl/gsl>
#include <memory>
#include <type_traits>

#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/sam.h>

#include <xoos/types/int.h>

namespace xoos::io {

struct SamFileDeleter {
  void operator()(samFile* h) const;
};

using SamFilePtr = std::unique_ptr<samFile, SamFileDeleter>;

struct SamHdrDeleter {
  void operator()(sam_hdr_t* h) const;
};

using SamHdrPtr = std::unique_ptr<sam_hdr_t, SamHdrDeleter>;

struct BamHdrDeleter {
  void operator()(bam_hdr_t* h) const;
};

using BamHdrPtr = std::unique_ptr<bam_hdr_t, BamHdrDeleter>;

struct Bam1Deleter {
  void operator()(bam1_t* b) const;
};

using Bam1Ptr = std::unique_ptr<bam1_t, Bam1Deleter>;

std::shared_ptr<bam1_t> MakeBam1SharedPtr(bam1_t* ptr);

using Bam1SharedPtr = std::shared_ptr<bam1_t>;

struct FastaIdxDeleter {
  void operator()(faidx_t* f) const;
};

using FaIdxPtr = std::unique_ptr<faidx_t, FastaIdxDeleter>;

struct HtsIdxDeleter {
  void operator()(hts_idx_t* i) const;
};

using HtsIdxPtr = std::unique_ptr<hts_idx_t, HtsIdxDeleter>;

struct HtsItrDeleter {
  void operator()(hts_itr_t* itr) const;
};

using HtsItrPtr = std::unique_ptr<hts_itr_t, HtsItrDeleter>;

struct HtsItrMultiDeleter {
  void operator()(hts_itr_multi_t* itr) const;
};

using HtsItrMultiPtr = std::unique_ptr<hts_itr_multi_t, HtsItrMultiDeleter>;

struct HtsFileDeleter {
  void operator()(htsFile* f) const;
};

using HtsFilePtr = std::unique_ptr<htsFile, HtsFileDeleter>;

class HtsRegListDeleter {
 public:
  explicit HtsRegListDeleter(s32 num_i);

  void operator()(gsl::owner<hts_reglist_t*> reglist) const;

 private:
  // Number of elements to free; tracked so the reglist is freed correctly.
  s32 _num_i{};
};

using HtsRegListPtr = std::unique_ptr<hts_reglist_t, HtsRegListDeleter>;

struct BamPlpDeleter {
  void operator()(bam_plp_t p) const;
};

using BamPlpPtr = std::unique_ptr<std::remove_pointer_t<bam_plp_t>, BamPlpDeleter>;

struct BgzfDeleter {
  void operator()(BGZF* f) const;
};

/// Owning smart pointer for a BGZF handle. Calls bgzf_close on destruction.
using BgzfPtr = std::unique_ptr<BGZF, BgzfDeleter>;

}  // namespace xoos::io
