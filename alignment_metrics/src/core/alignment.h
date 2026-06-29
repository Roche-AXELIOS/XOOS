#pragma once

#include <span>

#include <ankerl/unordered_dense.h>
#include <htslib/sam.h>

#include <xoos/types/int.h>

#include "alignment-metrics-options.h"
#include "core/hp-error-for-read.h"
#include "metadata/alignment-metadata.h"
#include "util/homopolymer-util.h"

namespace xoos::alignment_metrics {

/**
 * @brief Represents a CIGAR operation with its type and length.
 *
 * This struct is used to return CIGAR operation information from the Alignment class.
 * It provides a convenient way to access both the operation type (e.g., match, insertion, deletion)
 * and its length in a single structure, accounting for any trimming that may have been applied.
 */
struct Cigar {
  u32 type{};
  u32 length{};
};

/**
 * @brief Represents a read alignment with support for trimming and metric calculation.
 *
 * The Alignment class wraps a BAM alignment record and provides a view of the alignment
 * that can be modified (via trimming) without altering the underlying BAM record.
 * It stores parsed metadata, quality information, and CIGAR operations, and maintains
 * positions and lengths that can be adjusted when bases are trimmed from either end.
 *
 * This design allows efficient metric calculation on different views of the same alignment
 * without the overhead of creating new BAM records.
 */
class Alignment {
 public:
  // A non-owning pointer to the BAM alignment.
  bam1_t* alignment_ptr{};
  // Metadata for the alignment parsed from the read name and flags.
  AlignmentMetadata alignment_metadata{};
  // Base qualities and types for the read, including modified qualities/types based on homopolymer masking.
  ReadQualities qualities{};
  // Homopolymer error profile for the read, keyed by reference position of start position of reference homopolymers
  // covered by the read.
  ankerl::unordered_dense::map<s64, HomopolymerErrorProfileForRead> hp_error_profile{};

  /**
   * CIGAR operations for the alignment.
   * This is a non-owning span pointing to the CIGAR operations in the BAM alignment.
   * If the alignment is trimmed, this span will point to the trimmed CIGAR operations.
   *
   * NOTE: If an alignment is trimmed, the CIGAR operations should only be accessed via this member
   * and not directly from the BAM alignment pointer.
   */
  std::span<u32> cigars{};
  /**
   * The reference position of the alignment (0-based).
   * This is the position of the first aligned base on the reference.
   * If the alignment is trimmed, this position will be updated to reflect the new position after trimming.
   */
  s64 rpos{};
  /**
   * The query position of the alignment (0-based).
   * This is the position of the first aligned base on the read.
   * If the alignment is trimmed, this position will be updated to reflect the new position after trimming.
   */
  u32 qpos{};
  /**
   * The number of reference positions spanned by the alignment (including deletions, but not softclips or insertions).
   */
  u32 reference_length{};
  /**
   * Length of the read including softclips.
   * If the alignment is trimmed, this length will be updated to reflect the new length after trimming.
   * which should be equal to `alignment_ptr->core.l_qseq - number of trimmed bases`.
   */
  u32 read_length{};
  /**
   * Length of the read excluding softclips.
   */
  u32 read_length_without_softclips{};
  /**
   * Length of the first CIGAR operation.
   * This is needed because the first CIGAR operation may be partially trimmed.
   */
  u32 first_cigar_len{};
  /**
   * Length of the last CIGAR operation.
   * This is needed because the last CIGAR operation may be partially trimmed.
   */
  u32 last_cigar_len{};

  /**
   * Constructor to create an Alignment object from a BAM alignment pointer and options.
   * This will parse the alignment metadata and qualities, and initialize the CIGAR operations,
   * reference position, query position, reference span, and read lengths.
   */
  Alignment(bam1_t* ptr, const AlignmentMetricsOptions& options);

  /**
   * Get the length of the CIGAR operation at the given index.
   * This takes into account any trimming that may have occurred.
   * If an alignment has been trimmed, this function should be used instead of directly accessing the CIGAR operations.
   */
  u32 GetCigarLen(size_t cigar_index) const;

  /**
   * Get the CIGAR operation and its length at the given index.
   * This takes into account any trimming that may have occurred.
   * If an alignment has been trimmed, this function should be used instead of directly accessing the CIGAR operations.
   */
  Cigar GetCigarAt(size_t cigar_index) const;

  /**
   * Trim the alignment by removing the specified number of bases from the start and end of the read.
   * This will update the CIGAR operations, reference position, query position, reference span, and read lengths
   * to reflect the trimmed state.
   * An alignment can only be trimmed once. If an attempt is made to trim an already trimmed alignment,
   * an exception will be thrown.
   */
  void Trim(u32 trim_leading_bases, u32 trim_trailing_bases);

  u32 CountAlignedBases() const;

 private:
  /**
   * Indicates whether the alignment has been trimmed.
   * If true, the alignment has been trimmed and the CIGAR operations,
   * positions, and lengths reflect the trimmed state. An alignment with `trimmed` set to true
   * may have a different `rpos`, `reference_length`, `read_length`, and `read_length_without_softclips` compared to
   * the original alignment and cannot be trimmed again.
   *
   * If false, the alignment is in its original state.
   *
   * Private to prevent accidental modification.
   */
  bool _trimmed{};
};

}  // namespace xoos::alignment_metrics
