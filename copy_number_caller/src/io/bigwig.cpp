#include "io/bigwig.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <regex>
#include <stdexcept>
#include <string>
#include <string_view>

#include <libbigwig/bigWig.h>

#include <xoos/error/error.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/adopt.h>
#include <xoos/util/file-functions.h>

#include "io/fai.h"

namespace xoos::cnc {
void BigWigDeleter::operator()(bigWigFile_t* const fp) const {
  if (fp != nullptr) {
    bwClose(fp);
  }
}

void BwOverlappingIntervalsDeleter::operator()(bwOverlappingIntervals_t* const intervals) const {
  if (intervals != nullptr) {
    bwDestroyOverlappingIntervals(intervals);
  }
}

namespace {

constexpr std::string_view kReferenceMismatchMessage =
    "The mappability BigWig file is incompatible with the provided reference.";
constexpr size_t kBigWigStatsBinCount = 1;

// Mappability tracks may omit non-canonical contigs that are not used for interval generation.
const std::regex kNonCanonicalContigPattern("(_|chrEBV|chrM)");

bool IsCanonicalContig(const std::string& contig) {
  return !std::regex_search(contig, kNonCanonicalContigPattern);
}

BigWigPtr OpenBigWigFile(const fs::path& fname) {
  return AdoptUnique<bigWigFile_t, BigWigDeleter>(bwOpen(fname.c_str(), nullptr, "r"));  // NOSONAR(M23_286): C adopt
}

BigWig OpenMappabilityBigWig(const fs::path& fname) {
  try {
    return BigWig(fname);
  } catch (const std::runtime_error& e) {
    throw error::Error("Could not open mappability BigWig file {}: {}", fname.string(), e.what());
  }
}

u32 ToBigWigCoordinate(const size_t value, const std::string_view coordinate_name) {
  constexpr auto kMaxBigWigCoordinate = static_cast<size_t>(std::numeric_limits<u32>::max());
  if (value > kMaxBigWigCoordinate) {
    throw error::Error("BigWig {} coordinate {} exceeds the supported uint32 range", coordinate_name, value);
  }
  return static_cast<u32>(value);
}

}  // namespace

BigWig::BigWig(const fs::path& fname) : _fname(fname), _fp(OpenBigWigFile(fname)) {
  if (_fp == nullptr) {
    throw error::Error(
        "The provided BigWig (.bw) file is invalid due to corruption, truncation, "
        "or unsupported format. Please check the file's integrity and format.");
  }
}

BigWig::~BigWig() {
  bwCleanup();
}

/**
 * @brief get mean value over the entire interval
 * @param contig
 * @param start
 * @param end
 * @return
 */
f64 BigWig::GetMean(const std::string_view contig, const size_t start, const size_t end) const {
  const std::string contig_string(contig);
  const u32 start_u32 = ToBigWigCoordinate(start, "start");
  const u32 end_u32 = ToBigWigCoordinate(end, "end");
  const BwStatsPtr stats(  // NOSONAR(M23_286): C adopt
      bwStats(_fp.get(), contig_string.c_str(), start_u32, end_u32, kBigWigStatsBinCount, mean));
  if (!stats) {
    throw error::Error("Could not retrieve region {}:{}-{} from BigWig", contig, start, end);
  }
  return *stats;
}

BwOverlappingIntervalsPtr BigWig::GetInterval(const std::string_view contig,
                                              const size_t start,
                                              const size_t end) const {
  const std::string contig_string(contig);
  const u32 start_u32 = ToBigWigCoordinate(start, "start");
  const u32 end_u32 = ToBigWigCoordinate(end, "end");
  return AdoptUnique<bwOverlappingIntervals_t, BwOverlappingIntervalsDeleter>(
      bwGetOverlappingIntervals(_fp.get(), contig_string.c_str(), start_u32, end_u32));  // NOSONAR(M23_286): C adopt
}

vec<std::string> BigWig::GetHeaderContigs() const {
  if (_fp == nullptr || _fp->cl == nullptr) {
    throw error::Error("The provided BigWig (.bw) file is missing chromosome metadata: {}", _fname.string());
  }

  vec<std::string> header_contigs;
  header_contigs.reserve(static_cast<size_t>(_fp->cl->nKeys));
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic) - indexing libBigWig chromosome C arrays
  for (int64_t i = 0; i < _fp->cl->nKeys; ++i) {
    if (_fp->cl->chrom[i] == nullptr) {
      throw error::Error("The provided BigWig (.bw) file has invalid chromosome metadata: {}", _fname.string());
    }
    header_contigs.emplace_back(_fp->cl->chrom[i]);
  }
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  return header_contigs;
}

bool BigWig::CanQueryContig(const std::string_view contig) const {
  // A contig is queryable only if libBigWig can resolve its name AND retrieve summary data over the
  // whole contig. Missing name -> bwStats returns NULL; present name but unreachable data (blocks not
  // sorted to match the contig index) -> bwStats returns NaN. Both mean CNC's run-time lookups would
  // fail, so treat either as not queryable. The whole-contig query (rather than a small window) avoids
  // false negatives from valid gaps such as telomeric N runs, and exercises the same path as GetMean.
  const std::string contig_string(contig);
  const u32 tid = bwGetTid(_fp.get(), contig_string.c_str());
  if (tid == std::numeric_limits<u32>::max()) {
    return false;
  }
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic) - indexing libBigWig length C array by tid
  const u32 contig_length = _fp->cl->len[tid];
  const BwStatsPtr stats(  // NOSONAR(M23_286): C adopt
      bwStats(_fp.get(), contig_string.c_str(), 0, contig_length, kBigWigStatsBinCount, mean));
  return stats && std::isfinite(*stats);
}

void BigWig::ValidateReferenceContigs(const vec<std::string>& reference_contigs) const {
  // Touch header metadata up front so a structurally broken file fails with a clear message.
  static_cast<void>(GetHeaderContigs());

  for (const auto& contig : reference_contigs) {
    if (!IsCanonicalContig(contig)) {
      continue;
    }
    if (!CanQueryContig(contig)) {
      throw error::Error(
          "{} File: {}. The reference contig {} could not be queried in the mappability BigWig. The contig is either "
          "missing from the file, or its data is unreachable because the BigWig was not sorted to match its contig "
          "index.",
          kReferenceMismatchMessage,
          _fname.string(),
          contig);
    }
  }
}

BigWig ValidateAndOpenMappabilityBigWig(const fs::path& mappability_bigwig_fname,
                                        const fs::path& reference_genome_fai_fname) {
  file::CheckFileIsReadable(mappability_bigwig_fname);
  try {
    if (fs::file_size(mappability_bigwig_fname) == 0) {
      throw error::Error("Mappability file is empty: {}", mappability_bigwig_fname.string());
    }
  } catch (const fs::filesystem_error& e) {
    throw error::Error("Could not inspect mappability file {}: {}", mappability_bigwig_fname.string(), e.what());
  }

  // The mappability BigWig must be queryable for every canonical reference contig CNC uses. Contig
  // order in the BigWig header does not need to match the reference: libBigWig resolves contigs by
  // name, so ordering is validated implicitly by confirming each contig can be queried.
  const vec<std::string> reference_contigs = GetContigsInOrder(reference_genome_fai_fname);
  BigWig mapp_bigwig = OpenMappabilityBigWig(mappability_bigwig_fname);
  mapp_bigwig.ValidateReferenceContigs(reference_contigs);
  return mapp_bigwig;
}

}  // namespace xoos::cnc
