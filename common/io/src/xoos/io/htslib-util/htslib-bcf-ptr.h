/**
 * @file htslib-bcf-ptr.h
 * @brief RAII smart pointers for HTSlib BCF/VCF/tabix handles.
 *
 * Mirrors htslib-ptr.h (which covers SAM/BAM handles) for the BCF/VCF side.
 * Each deleter null-checks before calling the matching htslib destroy function.
 */
#pragma once

#include <memory>

#include <htslib/tbx.h>
#include <htslib/vcf.h>

namespace xoos::io {

struct BcfHdrDeleter {
  void operator()(bcf_hdr_t* h) const;
};

/// Owning smart pointer for a bcf_hdr_t. Calls bcf_hdr_destroy on destruction.
using BcfHdrPtr = std::unique_ptr<bcf_hdr_t, BcfHdrDeleter>;

struct Bcf1Deleter {
  void operator()(bcf1_t* b) const;
};

/// Owning smart pointer for a bcf1_t record. Calls bcf_destroy on destruction.
using Bcf1Ptr = std::unique_ptr<bcf1_t, Bcf1Deleter>;

struct TbxDeleter {
  void operator()(tbx_t* t) const;
};

/// Owning smart pointer for a tabix index. Calls tbx_destroy on destruction.
using TbxPtr = std::unique_ptr<tbx_t, TbxDeleter>;

}  // namespace xoos::io
