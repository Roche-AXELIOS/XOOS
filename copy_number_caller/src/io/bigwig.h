#pragma once

#include <cstddef>
#include <memory>

#include <xoos/types/float.h>
extern "C" {
#include <libbigwig/bigWig.h>
}

#include <string>
#include <string_view>

#include <xoos/io/malloc-ptr.h>
#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {

struct BigWigDeleter {
  void operator()(bigWigFile_t* fp) const;
};

using BigWigPtr = std::unique_ptr<bigWigFile_t, BigWigDeleter>;

using BwStatsPtr = io::MallocPtr<f64>;

struct BwOverlappingIntervalsDeleter {
  void operator()(bwOverlappingIntervals_t* intervals) const;
};

using BwOverlappingIntervalsPtr = std::unique_ptr<bwOverlappingIntervals_t, BwOverlappingIntervalsDeleter>;

class BigWig {
 public:
  BigWig() = delete;

  ~BigWig();

  BigWig(const BigWig&) = delete;
  BigWig(BigWig&&) noexcept = default;
  BigWig& operator=(const BigWig&) = delete;
  BigWig& operator=(BigWig&&) noexcept = default;

  explicit BigWig(const fs::path& fname);

  f64 GetMean(std::string_view contig, size_t start, size_t end) const;
  BwOverlappingIntervalsPtr GetInterval(std::string_view contig, size_t start, size_t end) const;

  /// @brief Confirm that every canonical reference contig can actually be queried in this BigWig.
  ///
  /// A file that opens cleanly can still fail later if a needed contig is absent, or if its data is
  /// unreachable because the source bedGraph was not sorted to match the BigWig contig index. This
  /// probes each canonical reference contig and throws if any is missing or has no retrievable data,
  /// turning a mid-run "could not retrieve region" into an upfront error naming the offending contig.
  ///
  /// Only canonical contigs are checked: alt, decoy, EBV, and mitochondrial contigs are routinely
  /// absent from a mappability track and excluded by CNC's interval generation, so requiring them
  /// would reject valid files.
  void ValidateReferenceContigs(const vec<std::string>& reference_contigs) const;

 private:
  vec<std::string> GetHeaderContigs() const;
  bool CanQueryContig(std::string_view contig) const;

  fs::path _fname;
  BigWigPtr _fp;
};

BigWig ValidateAndOpenMappabilityBigWig(const fs::path& mappability_bigwig_fname,
                                        const fs::path& reference_genome_fai_fname);

}  // namespace xoos::cnc
