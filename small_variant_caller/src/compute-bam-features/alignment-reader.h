#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/types/vec.h>
#include <xoos/util/simple-cache.h>

namespace xoos::svc {

using SamHdrSharedPtr = std::shared_ptr<sam_hdr_t>;
using HtsIdxSharedPtr = std::shared_ptr<hts_idx_t>;

/**
 * @brief Container for alignment file components including file pointer, header, and index
 */
struct AlignmentReader {
  io::HtsFilePtr fp;
  SamHdrSharedPtr hdr;
  HtsIdxSharedPtr idx;
};

/**
 * @brief Thread-safe cache for alignment file headers and indices
 *
 * This class exists to share read-only header and index data across multiple threads
 * to avoid loading the same data into memory multiple times. When multiple threads
 * need to access the same alignment file, they can share the header and index objects
 * while maintaining separate file pointers for thread safety.
 *
 * The cache uses weak references internally, so cached data is automatically cleaned
 * up when no longer referenced by any AlignmentReader instances.
 */
class AlignmentReaderCache {
 public:
  /**
   * @brief Opens a single alignment file and creates an AlignmentReader
   * @param fn Path to the alignment file
   * @param mode File open mode (e.g., "r" for read)
   * @return AlignmentReader containing file pointer, header, and index
   */
  AlignmentReader Open(const fs::path& fn, const std::string& mode);

  /**
   * @brief Opens multiple alignment files and creates a vector of AlignmentReaders
   * @param fn Vector of paths to alignment files
   * @param mode File open mode (e.g., "r" for read)
   * @return Vector of AlignmentReaders for all input files
   */
  vec<AlignmentReader> Open(const vec<fs::path>& fn, const std::string& mode);

 private:
  /**
   * @brief Extracts file path from htsFile pointer for cache key transformation
   * @param fp Pointer to htsFile
   * @return File system path extracted from the htsFile
   */
  static fs::path HtsFileToPath(htsFile* fp);

  /**
   * @brief Reads SAM header from file and wraps in shared_ptr
   * @param fp Pointer to htsFile to read header from
   * @return Shared pointer to loaded SAM header
   */
  static SamHdrSharedPtr SamHdrRead(htsFile* fp);

  /**
   * @brief Loads SAM index from file and wraps in shared_ptr
   * @param fp Pointer to htsFile to load index for
   * @return Shared pointer to loaded SAM index
   */
  static HtsIdxSharedPtr SamIndexLoad(htsFile* fp);

  // Cache for SAM headers keyed by file path to enable sharing across threads
  util::SimpleCache<htsFile*, SamHdrSharedPtr::element_type, fs::path> _hdr_cache{HtsFileToPath, SamHdrRead};

  // Cache for SAM indices keyed by file path to enable sharing across threads
  util::SimpleCache<htsFile*, HtsIdxSharedPtr::element_type, fs::path> _idx_cache{HtsFileToPath, SamIndexLoad};
};

}  // namespace xoos::svc
