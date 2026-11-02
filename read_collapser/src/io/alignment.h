#pragma once

#include <memory>

#include <htslib/sam.h>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/read-name/umi.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/yc-decode/yc-decoder.h>

namespace xoos::read_collapser {

struct Cluster;

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

/**
 * An Alignment represents a single read alignment, it tracks
 * the original record, BAM record as well as additional computed information.
 */
struct Alignment {
  /// The original record read from the BAM file, this is immutable and should not be modified.
  io::Bam1Ptr record;
  /// A non-owning pointer to the cluster this alignment belongs to.
  Cluster* cluster{};
  /// The UMI sequence at the 5' end of the read.
  read_name::Umi umi5p;
  /// The UMI sequence at the 3' end of the read.
  read_name::Umi umi3p;
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
  /**
   * Whether the read is complete and untruncated at the 3' end.
   */
  bool is_3p_complete{};
  /**
   * Whether the read is complete and untruncated at the 5' end.
   */
  bool is_5p_complete{};
  /**
   * Whether the read is a hairpin duplex read (e.g. SBX-D or SBX-FAST)
   */
  bool is_hairpin_read{};
  // The exclusive end position of the read on the reference genome. We cache this value to avoid repeatedly iterating
  // over the CIGAR to calculate it.
  u32 end_pos{};
  // The original flags of the read before any modifications during deconvolution.
  u32 original_flags{};

  // Create an Alignment object from a BAM record, parsing the read name to extract UMI and adapter information if
  // specified.
  explicit Alignment(io::Bam1Ptr record, bool expect_umis, bool ignore_read_name_parsing_errors);

  // Create an Alignment object from a BAM record, without parsing the read name for UMI and adapter information.
  explicit Alignment(io::Bam1Ptr record);

  // A read is partial if the insert is incomplete. A read can be truncated on the 5' end, 3' end, or both ends.
  // In all cases, the read is considered partial.
  bool IsPartial() const;

  /// Returns true if the read is complete and untruncated at the 5' end.
  bool IsFivePrimeComplete() const;

  /// Returns true if the read is complete and untruncated at the 3' end.
  bool IsThreePrimeComplete() const;

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

  bool HasSaTag() const;
};

using AlignmentPtr = std::shared_ptr<Alignment>;

/**
 * @brief Parses the UMI (Unique Molecular Identifier) and read adapter information from the read names of
 *        the given alignment and updates the alignment object with the parsed UMI values.
 *
 * @param alignment A reference to the alignment object to be processed.
 * @param expect_umis A boolean indicating whether UMIs are expected in the read names. If true, the function will
 *                    attempt to parse UMIs from the read names and update the alignment objects accordingly.
 * @param ignore_parsing_errors A boolean indicating whether to ignore parsing errors. If true, the function will
 *                              set UMI values to std::nullopt for alignments with parsing errors instead of throwing
 *                              exceptions. If false, the function will throw exceptions for any parsing errors
 *                              encountered.
 *
 * @throws error::Error If the UMI cannot be parsed from the read name and `ignore_parsing_errors` is false.
 */
void UpdateReadAdapterInfo(Alignment& alignment, bool expect_umis, bool ignore_parsing_errors);

}  // namespace xoos::read_collapser
