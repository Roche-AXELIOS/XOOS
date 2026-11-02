#include "core/pileup.h"

#include <algorithm>
#include <string>

#include <xoos/types/int.h>
#include <xoos/util/math.h>

#include "core/error-counts.h"
#include "core/hp-error-for-read.h"
#include "core/substitution-lookup.h"
#include "core/super-region.h"
#include "metadata/alignment-metadata.h"
#include "metrics/accuracy-metrics/hp-stat.h"
#include "metrics/metrics.h"

namespace xoos::alignment_metrics {

void DepthProfile::Reset() {
  _depth = 0;
  _concordant_duplex_depth = 0;
  _filtered_depth = 0;
  has_coverage = false;
}

void DepthProfile::AddDepth(const bool is_concordant_duplex, const bool passes_filters) {
  ++_depth;
  if (is_concordant_duplex) {
    ++_concordant_duplex_depth;
  }
  if (passes_filters) {
    ++_filtered_depth;
  }
}

u32 DepthProfile::GetDepth() const {
  return _depth;
}

u32 DepthProfile::GetConcordantDuplexDepth() const {
  return _concordant_duplex_depth;
}

u32 DepthProfile::GetFilteredDepth() const {
  return _filtered_depth;
}

void ErrorProfile::Reset() {
  _total_errors = CountsByErrorTypeAndReadProperty<u32>{};
  _errors_by_cluster_size.clear();
  _errors_by_substitution_type.clear();
  _alt_counts.clear();
  _qscore_stats.clear();
  _inserted_bases = 0;
  _deleted_bases = 0;
}

void ErrorProfile::AddAlignedBase(const AlignmentMetadata& alignment_metadata) {
  if (alignment_metadata.cluster_size > 0) {
    if (alignment_metadata.cluster_size >= _errors_by_cluster_size.size()) {
      _errors_by_cluster_size.resize(alignment_metadata.cluster_size + 1);
    }
    const size_t cluster_size_index = math::SatSub(alignment_metadata.cluster_size, u32{1});
    auto& cluster_stats = _errors_by_cluster_size[cluster_size_index];
    cluster_stats.totals.Add(alignment_metadata);
  }
  _total_errors.totals.Add(alignment_metadata);
}

void ErrorProfile::AddSubstitution(const char alt_base,
                                   const char ref_base,
                                   const AlignmentMetadata& alignment_metadata) {
  if (alignment_metadata.cluster_size > 0) {
    if (alignment_metadata.cluster_size >= _errors_by_cluster_size.size()) {
      _errors_by_cluster_size.resize(alignment_metadata.cluster_size + 1);
    }
    const size_t cluster_size_index = math::SatSub(alignment_metadata.cluster_size, u32{1});
    auto& cluster_stats = _errors_by_cluster_size[cluster_size_index];
    cluster_stats.substitutions.Add(alignment_metadata);
  }
  ++_alt_counts[std::string{alt_base}].substitutions;
  _total_errors.substitutions.Add(alignment_metadata);
  _errors_by_substitution_type[SubstitutionLookup::GetIndexForSubstitution(ref_base, alt_base)].Add(alignment_metadata);
}

void ErrorProfile::AddInsertion(const std::string& ins_seq, const AlignmentMetadata& alignment_metadata) {
  if (alignment_metadata.cluster_size > 0) {
    if (alignment_metadata.cluster_size >= _errors_by_cluster_size.size()) {
      _errors_by_cluster_size.resize(alignment_metadata.cluster_size + 1);
    }
    const size_t cluster_size_index = math::SatSub(alignment_metadata.cluster_size, u32{1});
    auto& cluster_stats = _errors_by_cluster_size[cluster_size_index];
    cluster_stats.insertions.Add(alignment_metadata);
  }
  ++_alt_counts[ins_seq].insertions;
  _total_errors.insertions.Add(alignment_metadata);
  _inserted_bases += static_cast<u32>(ins_seq.length());
}

void ErrorProfile::AddDeletion(const std::string& del_seq, const AlignmentMetadata& alignment_metadata) {
  if (alignment_metadata.cluster_size > 0) {
    if (alignment_metadata.cluster_size >= _errors_by_cluster_size.size()) {
      _errors_by_cluster_size.resize(alignment_metadata.cluster_size + 1);
    }
    const size_t cluster_size_index = math::SatSub(alignment_metadata.cluster_size, u32{1});
    auto& cluster_stats = _errors_by_cluster_size[cluster_size_index];
    cluster_stats.deletions.Add(alignment_metadata);
  }
  ++_alt_counts[del_seq].deletions;
  _total_errors.deletions.Add(alignment_metadata);
  _deleted_bases += static_cast<u32>(del_seq.length());
}

void ErrorProfile::AddQscoreData(u8 qscore, bool is_mismatch, const AlignmentMetadata& alignment_metadata) {
  _qscore_stats[qscore].Add(alignment_metadata, is_mismatch);
}

u32 ErrorProfile::GetDepth() const {
  CountsByErrorType<u32> stats;
  for (const auto& [_, counts_by_error_type] : _alt_counts) {
    stats.insertions += counts_by_error_type.insertions;
    stats.deletions += counts_by_error_type.deletions;
  }
  return stats.insertions + stats.deletions + _total_errors.totals.total;
}

f64 ErrorProfile::GetMaxAltAlleleFraction() const {
  const u32 depth = GetDepth();
  u32 max_alt_count = 0;
  for (const auto& [_, counts_by_error_type] : _alt_counts) {
    max_alt_count = std::max(max_alt_count, counts_by_error_type.insertions);
    max_alt_count = std::max(max_alt_count, counts_by_error_type.deletions);
    max_alt_count = std::max(max_alt_count, counts_by_error_type.substitutions);
  }
  return depth > 0 ? static_cast<f64>(max_alt_count) / static_cast<f64>(depth) : 0.0;
}

const CountsByErrorTypeAndReadProperty<u32>& ErrorProfile::GetTotalErrors() const {
  return _total_errors;
}

const vec<CountsByErrorTypeAndReadProperty<u32>>& ErrorProfile::GetErrorsByClusterSize() const {
  return _errors_by_cluster_size;
}

const ankerl::unordered_dense::map<u16, CountsByReadProperty<u32>>& ErrorProfile::GetErrorsBySubstitutionType() const {
  return _errors_by_substitution_type;
}

const ankerl::unordered_dense::map<u8, MismatchCounts<u32>>& ErrorProfile::GetQscoreStats() const {
  return _qscore_stats;
}

const ankerl::unordered_dense::map<std::string, CountsByErrorType<u32>>& ErrorProfile::GetAltAlleleCounts() const {
  return _alt_counts;
}

u32 ErrorProfile::GetInsertedBases() const {
  return _inserted_bases;
}

u32 ErrorProfile::GetDeletedBases() const {
  return _deleted_bases;
}

u32 HomopolymerErrorProfile::HpGetAltCountAtPosition(const vec<HpErrorByPosition>& hp_errors_by_position, u8 position) {
  u32 alt_depth = 0;
  for (const auto& hp_error_by_position : hp_errors_by_position) {
    alt_depth += ((hp_error_by_position & (1U << position)) != 0U) ? 1 : 0;
  }
  return alt_depth;
}

void HomopolymerErrorProfile::HpClearErrorsAtPosition(vec<HpErrorByPosition>& hp_errors_by_position, u8 position) {
  for (auto& hp_error_by_position : hp_errors_by_position) {
    hp_error_by_position &= ~(1U << position);
  }
}

u32 HomopolymerErrorProfile::HpGetNumberOfReadsWithErrors(const vec<HpErrorByPosition>& hp_errors_by_position) {
  u32 count = 0;
  for (const auto& hp_error_by_position : hp_errors_by_position) {
    count += hp_error_by_position != 0 ? 1 : 0;
  }
  return count;
}

void HomopolymerErrorProfile::Reset() {
  _insertions_by_position.clear();
  _deletions_by_position.clear();
  _indels_by_position.clear();
  _substitutions_by_position.clear();
  total_reads = 0;
  spanning_reads = 0;
  discordant_reads = 0;
  low_quality_reads = 0;
  effective_reads = 0;
}

void HomopolymerErrorProfile::AddErrorCounts(const HomopolymerErrorProfileForRead& hp_error_profile_for_read) {
  if (hp_error_profile_for_read.GetInsertionByPosition() != 0) {
    _insertions_by_position.emplace_back(hp_error_profile_for_read.GetInsertionByPosition());
  }
  if (hp_error_profile_for_read.GetDeletionByPosition() != 0) {
    _deletions_by_position.emplace_back(hp_error_profile_for_read.GetDeletionByPosition());
  }
  if (hp_error_profile_for_read.GetSubstitutionByPosition() != 0) {
    _substitutions_by_position.emplace_back(hp_error_profile_for_read.GetSubstitutionByPosition());
  }
  if (hp_error_profile_for_read.GetIndelByPosition() != 0) {
    _indels_by_position.emplace_back(hp_error_profile_for_read.GetIndelByPosition());
  }
}

void HomopolymerErrorProfile::FilterErrorsByAltAlleleFraction(u32 min_depth,
                                                              f64 max_alt_allele_fraction,
                                                              u8 hp_length) {
  if (effective_reads >= min_depth) {
    const auto effective_reads_f64 = static_cast<f64>(effective_reads);
    for (auto i = 0; i < hp_length; ++i) {
      const f64 insertion_fraction_at_i = HpGetAltCountAtPosition(_insertions_by_position, i) / effective_reads_f64;
      if (insertion_fraction_at_i > max_alt_allele_fraction) {
        HpClearErrorsAtPosition(_insertions_by_position, i);
      }
      const f64 deletion_fraction_at_i = HpGetAltCountAtPosition(_deletions_by_position, i) / effective_reads_f64;
      if (deletion_fraction_at_i > max_alt_allele_fraction) {
        HpClearErrorsAtPosition(_deletions_by_position, i);
      }
      const f64 indel_fraction_at_i = HpGetAltCountAtPosition(_indels_by_position, i) / effective_reads_f64;
      if (indel_fraction_at_i > max_alt_allele_fraction) {
        HpClearErrorsAtPosition(_indels_by_position, i);
      }
      const f64 substitution_fraction_at_i =
          HpGetAltCountAtPosition(_substitutions_by_position, i) / effective_reads_f64;
      if (substitution_fraction_at_i > max_alt_allele_fraction) {
        HpClearErrorsAtPosition(_substitutions_by_position, i);
      }
    }
  }
}

u32 HomopolymerErrorProfile::GetSubstitutionCount() const {
  return HpGetNumberOfReadsWithErrors(_substitutions_by_position);
}

u32 HomopolymerErrorProfile::GetInsertionCount() const {
  return HpGetNumberOfReadsWithErrors(_insertions_by_position);
}

u32 HomopolymerErrorProfile::GetDeletionCount() const {
  return HpGetNumberOfReadsWithErrors(_deletions_by_position);
}

u32 HomopolymerErrorProfile::GetIndelCount() const {
  return HpGetNumberOfReadsWithErrors(_indels_by_position);
}

void ReferenceContext::Reset() {
  hp_base = kHpStatUninitializedBase;
  hp_length = 0;
  relative_position_within_hp = 0;
  window_gc_percentage = -1;
}

namespace {

// Compute the mean-coverage histogram bin for a subregion from its per-position filtered depths.
u32 SubregionMeanCoverageBin(const vec<Pileup>& pileups, const SuperRegion& super_region, const Interval& region) {
  u64 total_aligned_bases = 0;
  const u64 region_length = ToUnsigned(region.end - region.start);
  for (s64 pos = region.start; pos < region.end; ++pos) {
    const auto pileup_index = static_cast<size_t>(pos - super_region.start);
    if (pileup_index < pileups.size() && pileups.at(pileup_index).valid) {
      total_aligned_bases += pileups.at(pileup_index).depth_profile.GetFilteredDepth();
    }
  }
  const f64 mean_coverage =
      region_length > 0 ? static_cast<f64>(total_aligned_bases) / static_cast<f64>(region_length) : 0.0;
  return static_cast<u32>(std::floor(mean_coverage));
}

// Record a subregion's mean-coverage bin into a single Metrics container's uniformity histogram.
void RecordMeanCoverageBin(Metrics& metrics, const u32 mean_coverage_bin) {
  if (metrics.coverage_metrics.has_value() && metrics.coverage_metrics->coverage_uniformity_metrics.has_value()) {
    ++metrics.coverage_metrics->coverage_uniformity_metrics->region_counts_by_mean_coverage[mean_coverage_bin];
  }
}

// Adds a valid pileup's depth into the coverage histograms.
void AggregateCoverage(const Pileup& pileup, Metrics& metrics) {
  if (pileup.depth_profile.has_coverage) {
    metrics.coverage_metrics->AddHistogramData(pileup.depth_profile.GetDepth(),
                                               pileup.depth_profile.GetConcordantDuplexDepth(),
                                               pileup.depth_profile.GetFilteredDepth(),
                                               pileup.is_part_of_hp);
  }
}

// Adds a valid pileup's error profile into the accuracy metrics, gated by depth and alt-allele-fraction thresholds.
void AggregateAccuracy(const Pileup& pileup, Metrics& metrics, const u32 min_depth, const f64 max_alt_allele_fraction) {
  auto& accuracy = metrics.accuracy_metrics.value();
  if (pileup.error_profile.GetDepth() < min_depth) {
    ++accuracy.base_level_accuracy_summary.positions_skipped_due_to_low_depth;
    ++accuracy.base_level_accuracy_summary.total_positions;
    return;
  }
  if (pileup.error_profile.GetMaxAltAlleleFraction() > max_alt_allele_fraction) {
    ++accuracy.base_level_accuracy_summary.positions_skipped_due_to_high_alt_allele_fraction;
    ++accuracy.base_level_accuracy_summary.total_positions;
    return;
  }
  accuracy.base_level_accuracy_summary += pileup.error_profile.GetTotalErrors();
  accuracy.base_level_accuracy_summary.inserted_bases += pileup.error_profile.GetInsertedBases();
  accuracy.base_level_accuracy_summary.deleted_bases += pileup.error_profile.GetDeletedBases();
  accuracy.errors_by_substitution_type += pileup.error_profile.GetErrorsBySubstitutionType();
  accuracy.qscore_stats += pileup.error_profile.GetQscoreStats();
  if (accuracy.errors_by_read_type.has_value()) {
    accuracy.errors_by_read_type.value() += pileup.error_profile.GetTotalErrors();
  }
  if (accuracy.errors_by_cluster_size.has_value()) {
    accuracy.errors_by_cluster_size.value() += pileup.error_profile.GetErrorsByClusterSize();
  }
  ++accuracy.base_level_accuracy_summary.total_positions;
}

// Adds a valid homopolymer-anchor pileup's homopolymer error profile into the homopolymer accuracy stats.
void AggregateHomopolymerStats(const Pileup& pileup,
                               Metrics& metrics,
                               const u32 min_depth,
                               const f64 max_alt_allele_fraction) {
  if (!pileup.is_part_of_hp || pileup.reference_context.relative_position_within_hp != 0) {
    return;
  }
  auto hp_error_profile = pileup.hp_error_profile;
  auto& hp_metrics = metrics.accuracy_metrics->hp_stats->hp_stats[HpStatKey{
      .base = pileup.reference_context.hp_base, .hp_length = pileup.reference_context.hp_length}];
  hp_metrics.total_reads += hp_error_profile.total_reads;
  hp_metrics.spanning_reads += hp_error_profile.spanning_reads;
  hp_metrics.discordant_reads += hp_error_profile.discordant_reads;
  hp_metrics.low_quality_reads += hp_error_profile.low_quality_reads;
  hp_metrics.effective_reads += hp_error_profile.effective_reads;
  // Filter out positions with an alt fraction greater than the max alt allele fraction if the effective
  // coverage is greater than the minimum depth
  hp_error_profile.FilterErrorsByAltAlleleFraction(
      min_depth, max_alt_allele_fraction, pileup.reference_context.hp_length);
  // Instead of counting the total number of errors, we count the number of reads that have an insertion,
  // deletion, or substitution at some position within the homopolymer region.
  hp_metrics.count_ins += hp_error_profile.GetInsertionCount();
  hp_metrics.count_del += hp_error_profile.GetDeletionCount();
  hp_metrics.count_sub += hp_error_profile.GetSubstitutionCount();
  hp_metrics.count_indel += hp_error_profile.GetIndelCount();
}

// Aggregate a single valid pileup into one Metrics container.
void AggregatePileupInto(const Pileup& pileup,
                         Metrics& target,
                         const u32 min_depth,
                         const f64 max_alt_allele_fraction) {
  if (target.coverage_metrics.has_value()) {
    AggregateCoverage(pileup, target);
  }
  if (target.accuracy_metrics.has_value()) {
    AggregateAccuracy(pileup, target, min_depth, max_alt_allele_fraction);
    AggregateHomopolymerStats(pileup, target, min_depth, max_alt_allele_fraction);
  }
}

}  // namespace

void UpdateCoverageUniformityHistogram(const vec<Pileup>& pileups,
                                       const SuperRegion& super_region,
                                       StratifiedMetrics& metrics) {
  for (const auto& region : super_region.subregions) {
    const u32 mean_coverage_bin = SubregionMeanCoverageBin(pileups, super_region, region);
    // Route this subregion's mean-coverage bin into every label container it belongs to.
    metrics.ForEachContainer(
        region.labels, [mean_coverage_bin](Metrics& target) { RecordMeanCoverageBin(target, mean_coverage_bin); });
  }
}

void UpdateCoverageUniformityHistogram(const vec<Pileup>& pileups, const SuperRegion& super_region, Metrics& metrics) {
  for (const auto& region : super_region.subregions) {
    RecordMeanCoverageBin(metrics, SubregionMeanCoverageBin(pileups, super_region, region));
  }
}

void AggregateMetricsFromPileups(const vec<Pileup>& pileups,
                                 const SuperRegion& super_region,
                                 StratifiedMetrics& metrics,
                                 const u32 min_depth,
                                 const f64 max_alt_allele_fraction) {
  // Subregions are disjoint (enforced by CreateSuperRegion), so iterating their positions visits
  // each valid pileup at most once. Routing per subregion lets us read the label bit vector from the
  // subregion's Interval instead of storing it on every pileup position.
  for (const auto& subregion : super_region.subregions) {
    metrics.ForEachContainer(
        subregion.labels, [&subregion, &super_region, &pileups, min_depth, max_alt_allele_fraction](Metrics& target) {
          for (s64 pos = subregion.start; pos < subregion.end; ++pos) {
            const auto pileup_index = static_cast<size_t>(pos - super_region.start);
            if (pileup_index >= pileups.size()) {
              continue;
            }
            const auto& pileup = pileups.at(pileup_index);
            if (pileup.valid) {
              AggregatePileupInto(pileup, target, min_depth, max_alt_allele_fraction);
            }
          }
        });
  }
}

void AggregateMetricsFromPileups(const vec<Pileup>& pileups,
                                 Metrics& metrics,
                                 const u32 min_depth,
                                 const f64 max_alt_allele_fraction) {
  for (const auto& pileup : pileups) {
    if (!pileup.valid) {
      continue;
    }
    AggregatePileupInto(pileup, metrics, min_depth, max_alt_allele_fraction);
  }
}

}  // namespace xoos::alignment_metrics
