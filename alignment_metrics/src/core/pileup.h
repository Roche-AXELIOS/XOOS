#pragma once

#include <string>

#include <ankerl/unordered_dense.h>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/error-counts.h"
#include "core/hp-error-for-read.h"
#include "metadata/alignment-metadata.h"
#include "metrics/metrics.h"

namespace xoos::alignment_metrics {

struct Pileup;

/**
 * A class to track depth information at a specific reference position.
 */
class DepthProfile {
 public:
  /**
   * Reset all depth counts to their initial state.
   */
  void Reset();

  /**
   * Indicates whether the position associated with this depth profile should be considered for coverage calculations.
   */
  bool has_coverage{false};
  /**
   * Adds depth information for a base aligned at this position.
   * @param is_concordant_duplex Whether the base is from a concordant duplex read.
   * @param passes_filters Whether the base passes quality and type filters.
   */
  void AddDepth(bool is_concordant_duplex, bool passes_filters);

  u32 GetDepth() const;

  u32 GetConcordantDuplexDepth() const;

  u32 GetFilteredDepth() const;

  friend void AggregateMetricsFromPileups(const vec<Pileup>& pileups,
                                          Metrics& metrics,
                                          u32 min_depth,
                                          f64 max_alt_allele_fraction);
  friend void UpdateCoverageUniformityHistogram(const vec<Pileup>& pileups,
                                                const SuperRegion& super_region,
                                                Metrics& metrics);

 private:
  u32 _depth{};
  u32 _concordant_duplex_depth{};
  u32 _filtered_depth{};
};

/**
 * A class to track error information at a specific reference position, including substitutions and indels.
 */
class ErrorProfile {
  CountsByErrorTypeAndReadProperty<u32> _total_errors{};
  vec<CountsByErrorTypeAndReadProperty<u32>> _errors_by_cluster_size{};
  ankerl::unordered_dense::map<u16, CountsByReadProperty<u32>> _errors_by_substitution_type{};
  ankerl::unordered_dense::map<std::string, CountsByErrorType<u32>> _alt_counts{};
  ankerl::unordered_dense::map<u8, MismatchCounts<u32>> _qscore_stats{};
  u32 _inserted_bases{};
  u32 _deleted_bases{};

 public:
  /**
   * Reset all error counts and statistics to their initial state.
   */
  void Reset();

  /**
   * Add an aligned base to the error profile, updating total and cluster-specific stats.
   * @param alignment_metadata Metadata for the alignment containing the base.
   */
  void AddAlignedBase(const AlignmentMetadata& alignment_metadata);
  /**
   * Add a substitution error (reference mismatch) to the error profile, updating total and cluster-specific stats.
   * @param alt_base The observed base in the read.
   * @param ref_base The reference base.
   * @param alignment_metadata Metadata for the alignment containing the substitution.
   */
  void AddSubstitution(char alt_base, char ref_base, const AlignmentMetadata& alignment_metadata);
  /**
   * Add an insertion error to the error profile, updating total and cluster-specific stats.
   * @param ins_seq The inserted sequence.
   * @param alignment_metadata Metadata for the alignment containing the insertion.
   */
  void AddInsertion(const std::string& ins_seq, const AlignmentMetadata& alignment_metadata);
  /**
   * Add a deletion error to the error profile, updating total and cluster-specific stats.
   * @param del_seq The deleted sequence.
   * @param alignment_metadata Metadata for the alignment containing the deletion.
   */
  void AddDeletion(const std::string& del_seq, const AlignmentMetadata& alignment_metadata);
  /**
   * Add quality score data to the error profile.
   * @param qscore The quality score of the base.
   * @param is_mismatch Whether the base is a mismatch.
   * @param alignment_metadata Metadata for the alignment containing the base.
   */
  void AddQscoreData(u8 qscore, bool is_mismatch, const AlignmentMetadata& alignment_metadata);
  /**
   * Get the total depth (number of bases including indels) at this position.
   * @return The total depth.
   */
  u32 GetDepth() const;
  /**
   * Get the maximum alt (variant) allele fraction at the position associated with this error profile.
   */
  f64 GetMaxAltAlleleFraction() const;

  const CountsByErrorTypeAndReadProperty<u32>& GetTotalErrors() const;

  const ankerl::unordered_dense::map<u16, CountsByReadProperty<u32>>& GetErrorsBySubstitutionType() const;

  const ankerl::unordered_dense::map<u8, MismatchCounts<u32>>& GetQscoreStats() const;

  const ankerl::unordered_dense::map<std::string, CountsByErrorType<u32>>& GetAltAlleleCounts() const;

  u32 GetInsertedBases() const;

  u32 GetDeletedBases() const;

  friend void AggregateMetricsFromPileups(const vec<Pileup>& pileups,
                                          Metrics& metrics,
                                          u32 min_depth,
                                          f64 max_alt_allele_fraction);
};

/**
 * A class to track homopolymer error profiles at a specific reference position.
 * If a position is the first base of a homopolymer, an instance of this class is created to track errors
 * within the homopolymer region. For positions that are part of a homopolymer but not the first base, an instance
 * is still created but only the `relative_position_within_hp` field is set to allow back-referencing to the
 * first base of the homopolymer.
 */
class HomopolymerErrorProfile {
  vec<HpErrorByPosition> _insertions_by_position;
  vec<HpErrorByPosition> _deletions_by_position;
  vec<HpErrorByPosition> _substitutions_by_position;
  vec<HpErrorByPosition> _indels_by_position;

  /**
   * Get the number of reads with an alt allele (mismatch or indel with respect to the reference) at a specific position
   * within the homopolymer.
   * @param hp_errors_by_position Vector of bit vectors indicating error positions for each read.
   * @param position The position within the homopolymer to check for alt alleles.
   */
  static u32 HpGetAltCountAtPosition(const vec<HpErrorByPosition>& hp_errors_by_position, u8 position);
  /**
   * Clear error flags at a specific position within the homopolymer across all reads.
   * @param hp_errors_by_position Vector of bit vectors indicating error positions for each read.
   * @param position The position within the homopolymer to clear errors for.
   */
  static void HpClearErrorsAtPosition(vec<HpErrorByPosition>& hp_errors_by_position, u8 position);
  /**
   * Get the number of reads that have any type of error (insertion, deletion, or substitution) at any position within
   * the homopolymer.
   * @param hp_errors_by_position Vector of bit vectors indicating error positions for each read.
   */
  static u32 HpGetNumberOfReadsWithErrors(const vec<HpErrorByPosition>& hp_errors_by_position);

 public:
  /**
   * Reset all error counts and statistics to their initial state.
   */
  void Reset();

  /**
   * The base character of the homopolymer (e.g., 'A', 'C', 'G', 'T').
   * This should be initialized when a `HomopolymerErrorProfile` object is created.
   */
  char base{kHpStatUninitializedBase};
  /**
   * The length of the homopolymer.
   * This should be initialized when a `HomopolymerErrorProfile` object is created.
   */
  u8 hp_length{0};
  /**
   * The position of this profile within the homopolymer (0-indexed).
   */
  u8 relative_position_within_hp{0};
  /**
   * The number of reads that cover this homopolymer position.
   * This should be updated as reads are processed.
   */
  u32 total_reads{0};
  /**
   * The number of reads that span the entire homopolymer with sufficient anchor bases on both sides.
   * This should be updated as reads are processed.
   */
  u32 spanning_reads{0};
  /**
   * The number of spanning reads that cover the homopolymer region associated with this profile but have discordant
   * bases within the homopolymer. This should be updated when reads are processed.
   */
  u32 discordant_reads{0};
  /**
   * The number of spanning reads that cover the homopolymer region associated with this profile but have bases that do
   * not pass quality or type filters. This should be updated when reads are processed.
   */
  u32 low_quality_reads{0};
  /**
   * The number of spanning reads that cover the homopolymer region associated with this profile and have bases that
   * pass quality and type filters and do not have discordant bases within the homopolymer.
   */
  u32 effective_reads{0};
  /**
   * Update the error counts associated with this homopolymer based on the per-read error profile.
   * @param hp_error_profile_for_read The homopolymer error profile for a single read.
   */
  void AddErrorCounts(const HomopolymerErrorProfileForRead& hp_error_profile_for_read);
  /**
   * Filter out errors at positions within the homopolymer that have an alt allele fraction greater than the specified
   * threshold, provided the effective coverage meets the minimum depth requirement.
   * @param min_depth Minimum effective coverage required to apply filtering.
   * @param max_alt_allele_fraction Maximum alt allele fraction allowed to retain errors at a position.
   */
  void FilterErrorsByAltAlleleFraction(u32 min_depth, f64 max_alt_allele_fraction);
  /**
   * Get the number of reads with a substitution at any position within the homopolymer.
   * @return The number of reads with a substitution.
   */
  u32 GetSubstitutionCount() const;
  /**
   * Get the number of reads with an insertion at any position within the homopolymer.
   * @return The number of reads with an insertion.
   */
  u32 GetInsertionCount() const;
  /**
   * Get the number of reads with a deletion at any position within the homopolymer.
   * @return The number of reads with a deletion.
   */
  u32 GetDeletionCount() const;
  /**
   * Get the number of reads with an indel (insertion or deletion) at any position within the homopolymer.
   * @return The number of reads with an indel.
   */
  u32 GetIndelCount() const;

  friend void AggregateMetricsFromPileups(const vec<Pileup>& pileups,
                                          Metrics& metrics,
                                          u32 min_depth,
                                          f64 max_alt_allele_fraction);
};

/**
 * Represents a pileup at a specific reference position, including depth, error, and homopolymer error profiles.
 */
struct Pileup {
  /**
   * Indicates whether this pileup position should be considered for metric calculations.
   * This should be initialized when a vector of `Pileup` objects is created for a genomic region, and
   * positions that are part of the region of interest should be marked as valid.
   */
  bool valid{false};
  /**
   * Indicates whether the position associated with this pileup is part of a homopolymer region.
   * This should be initialized when a vector of `Pileup` objects is created for a genomic region, and
   * positions that are part of homopolymers should be marked accordingly.
   * NOTE: Although homopolymer metrics are only collected at the first base of a homopolymer,
   * this flag is set for all bases within the homopolymer.
   */
  bool is_part_of_hp{false};
  /**
   * The depth profile for this pileup position, tracking coverage information.
   */
  DepthProfile depth_profile;
  /**
   * The error profile for this pileup position, tracking substitution and indel error counts
   * and additional error statistics like error counts by cluster size and by substitution type.
   */
  ErrorProfile error_profile;
  /**
   * The homopolymer error profile for this pileup position, if the position is the first base of a homopolymer.
   * For positions that are part of a homopolymer but not the first base, all fields except
   * `relative_position_within_hp` are left uninitialized. The `relative_position_within_hp` field is set so that it can
   * be used as a back-reference to the first base of the homopolymer when processing reads.
   */
  HomopolymerErrorProfile hp_error_profile;
};

/**
 * Update the mean coverage histogram in the provided `Metrics` object based on the coverage data
 * from the given vector of `Pileup` objects and the associated `SuperRegion`. This histogram has
 * bins representing mean coverage levels and the values for each bin represent the number of
 * target regions (subregions) that fall within the corresponding mean coverage level.
 *
 * @param pileups The vector of `Pileup` objects containing coverage data.
 * @param super_region The `SuperRegion` object representing the genomic region associated with the pileups. Used to
 * identify target regions for histogram calculation.
 * @param metrics The `Metrics` object to update with mean coverage histogram data.
 */
void UpdateCoverageUniformityHistogram(const vec<Pileup>& pileups, const SuperRegion& super_region, Metrics& metrics);

/**
 * Aggregate metrics from a vector of pileups into the provided `Metrics` object.
 * @param pileups The vector of `Pileup` objects to aggregate metrics from.
 * @param metrics The `Metrics` object to update with aggregated data.
 * @param min_depth Minimum depth required to include a position in accuracy metrics aggregation.
 * @param max_alt_allele_fraction Maximum alt allele fraction allowed to include a position in accuracy metrics
 * aggregation.
 */
void AggregateMetricsFromPileups(const vec<Pileup>& pileups,
                                 Metrics& metrics,
                                 u32 min_depth,
                                 f64 max_alt_allele_fraction);

}  // namespace xoos::alignment_metrics
