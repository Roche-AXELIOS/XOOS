#pragma once

#include <memory>
#include <optional>
#include <string>

#include <htslib/sam.h>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/yc-decode/yc-decoder.h>

namespace xoos::read_collapser {

struct Cluster;

using Umi = std::optional<std::string>;

/**
 * A duplex read consists of two reads, R1 and R2, linked with a loop adapter,
 * forming a hairpin structure. Demux produces a single intramolecular consensus
 * read, which can be deconvolved back into the constituent R1 and R2 reads. The
 * duplex strand property indicates whether the read is R1, R2, or simplex (not deconvolved).
 */
enum class DuplexStrand : u8 {
  kR1,
  kR2,
  kSimplex
};

enum class ReadStrand : u8 {
  kFwd,
  kRev
};

const std::string kNoUmi = "*";  /// Represents no UMI in the read name.

/**
 * An Alignment represents a single read alignment, it tracks
 * the original record, BAM record as well as additional computed information.
 */
struct Alignment {
  /// The original record read from the BAM file, this is immutable and should not be modified.
  io::Bam1Ptr record{};
  /// A non-owning pointer to the cluster this alignment belongs to.
  Cluster* cluster{};
  /// The UMI sequence at the 5' end of the read.
  Umi umi5p{};
  /// The UMI sequence at the 3' end of the read.
  Umi umi3p{};
  /// The duplex strand of the read, if applicable.
  DuplexStrand duplex_strand{DuplexStrand::kSimplex};
  /// Whether this alignment is a duplicate, only set during duplicate marking.
  bool is_duplicate{};
  /**
   * Whether the mate of this alignment in a duplex pair has a left simplex overhang (only relevant
   * for deconvolved duplex reads and should always be `false` for simplex reads).
   * This is used to determine whether or not to pad the leading gaps with 'P' in the consensus matrix.
   *
   * For example, given a duplex read AAAAACCCCCCTCCCCCCC with YC tag 5-6E7-0,
   * R1 has a left overhang of 5 bases (AAAAA) so `mate_has_left_overhang` would be
   * set to true for R2 (since R1 is the mate of R2 in the same duplex pair).
   */
  bool mate_has_left_overhang{};
  /**
   * Whether the mate of this alignment in a duplex pair has a right simplex overhang (only relevant
   * for deconvolved duplex reads and should always be `false` for simplex reads).
   * This is used to determine whether or not to pad the trailing gaps with 'P' in the consensus matrix.
   *
   * For example, given a duplex read AAAAACCCCCCTCCCCCCC with YC tag 0-6E7-5,
   * R1 has a right overhang of 5 bases (CCCCC) so `mate_has_right_overhang` would be
   * set to true for R2 (since R1 is the mate of R2 in the same duplex pair).
   */
  bool mate_has_right_overhang{};
  // The exclusive end position of the read on the reference genome. We cache this value to avoid repeatedly iterating
  // over the CIGAR to calculate it.
  u32 end_pos{};
  // The original flags of the read before any modifications during deconvolution.
  u32 original_flags{};

  explicit Alignment(io::Bam1Ptr record);

  // A read is partial if the insert is complete. Currently we have no easy mechanism to determine this in RC; it is
  // additional metadata to be added upstream in future. For now we assume that if the data has UMI and one UMI is
  // missing then it is partial (for example v9.2 and SBX-S), and if there are no UMI at all then the read is full (for
  // example duplex).
  bool IsPartial() const;

  /// Returns true if the read is mapped to the reverse strand.
  bool IsReverse() const;

  /// Returns true if the read is mapped to the forward strand.
  bool IsForward() const;

  /// The inclusive start position of the read on the reference genome.
  u32 StartPos() const;

  /// The exclusive end position of the read on the reference genome.
  u32 EndPos() const;

  /**
   * Returns a vector with the same size as the read length, where each element is the base type
   * at the corresponding position in the read. If the read does not have a YC tag, this function returns an empty
   * vector.
   */
  vec<yc_decode::BaseType> GetBaseTypes() const;
};

using AlignmentPtr = std::shared_ptr<Alignment>;

}  // namespace xoos::read_collapser
