#include "alignment-reader.h"

namespace xoos::svc {

AlignmentReader AlignmentReaderCache::Open(const fs::path& fn, const std::string& mode) {
  // Open new file pointer (each thread gets its own for thread safety)
  auto fp = io::HtsOpen(fn, mode);

  // Get cached header and index (shared across threads for memory efficiency)
  auto hdr = _hdr_cache.Get(fp.get());
  auto idx = _idx_cache.Get(fp.get());

  return {std::move(fp), std::move(hdr), std::move(idx)};
}

vec<AlignmentReader> AlignmentReaderCache::Open(const vec<fs::path>& fn, const std::string& mode) {
  vec<AlignmentReader> alignment_readers;
  alignment_readers.reserve(fn.size());

  // Open each file, potentially sharing headers and indices from cache
  for (const auto& f : fn) {
    alignment_readers.emplace_back(Open(f, mode));
  }
  return alignment_readers;
}

fs::path AlignmentReaderCache::HtsFileToPath(htsFile* fp) {
  return fs::path{fp->fn};
}

SamHdrSharedPtr AlignmentReaderCache::SamHdrRead(htsFile* fp) {
  auto hdr = io::SamHdrRead(fp);
  return {std::move(hdr)};
}

HtsIdxSharedPtr AlignmentReaderCache::SamIndexLoad(htsFile* fp) {
  auto idx = io::SamIndexLoad(fp, fp->fn);
  return {std::move(idx)};
}

}  // namespace xoos::svc
