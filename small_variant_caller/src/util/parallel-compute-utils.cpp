#include "parallel-compute-utils.h"

#include <algorithm>
#include <limits>
#include <mutex>
#include <string_view>
#include <unordered_set>

#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/for_each.hpp>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/util/math.h>

#include "compute-bam-features/region.h"
#include "compute-vcf-features/compute-vcf-features.h"
#include "compute-vcf-features/vcf-header-util.h"
#include "core/vcf-fields.h"
#include "vcf-to-bed/vcf-to-bed.h"
#include "vcf-util.h"

namespace xoos::svc {

WorkerContext::WorkerContext(const std::string& vcf_file,
                             const AnnotationVcfPaths& annotation_vcfs,
                             const std::vector<fs::path>& bam_inputs,
                             AlignmentReaderCache& alignment_reader_cache)
    // Set up pointers to BAM readers, SAM headers, and BAM indexes for all input BAM files
    : filter_variants_reader{vcf_file},
      vcf_feature_extraction_reader{vcf_file},
      alignment_readers{alignment_reader_cache.Open(bam_inputs, "r")} {
  // Set up the population allele frequency VCF reader if provided
  if (annotation_vcfs.popaf.has_value()) {
    popaf_reader = io::VcfReader(annotation_vcfs.popaf.value());
  }
  // Set up the Panel of Normals VCF reader if provided
  if (annotation_vcfs.pon.has_value()) {
    pon_reader = io::VcfReader(annotation_vcfs.pon.value());
  }
  // Set up the catalog VCF reader if provided
  if (annotation_vcfs.catalog.has_value()) {
    catalog_reader = io::VcfReader(annotation_vcfs.catalog.value());
  }
}

/**
 * @brief Partition a VCF file into regions based on the number of threads available.
 * Each region will contain approximately the same number of variants, allowing for parallel processing.
 * @param vcf_file Path to the VCF file to partition.
 * @param threads Number of threads to use for parallel processing.
 * @return A vector of TargetRegion objects representing the partitioned regions.
 */
vec<TargetRegion> PartitionVcfRegions(const fs::path& vcf_file,
                                      const size_t threads,
                                      std::optional<ChromIntervalsMap> bed_regions) {
  vec<TargetRegion> partitioned_regions;
  const auto safe_threads = std::max<size_t>(threads, 1);
  const u32 left_pad = 1;
  const u32 right_pad = 1;
  ChromIntervalsMap target_regions;
  if (bed_regions.has_value()) {
    target_regions = bed_regions.value();
  }

  // Only collapse intervals if they are immediately adjacent to each other. So, each chromosome would have
  // approximately one interval for each variant. The only exception is large deletions, which have large intervals.
  const u32 collapse_dist = 0;

  // Use a single thread to extract variant positions. If the VCF file has many chromosomes, then the parallelized
  // version (1 parallel task per chromomosome) may cause this step to become very slow.
  auto chrom_to_intervals =
      ExtractVariantIntervalsSingleThreaded(vcf_file, target_regions, left_pad, right_pad, collapse_dist);
  // The goal here is to assign roughly same number of variants per parallel task in each chromosome. Since features
  // have already been computed, there is no need to create parallel tasks for small equal-sized regions across the
  // entire genome. Otherwise, many tasks would be created for regions that do not overlap any variants at all.
  for (auto& [chrom, intervals] : chrom_to_intervals) {
    if (!intervals.empty()) {
      const auto num_intervals = intervals.size();
      if (num_intervals > safe_threads) {
        // All partitions of this chromosome would contain near-equal number of intervals.
        const auto num_intervals_per_thread = (num_intervals + safe_threads - 1) / safe_threads;
        const auto last_i = num_intervals - 1;
        for (size_t i = 0; i < num_intervals; i += num_intervals_per_thread) {
          const auto max_i = std::min(last_i, i + num_intervals_per_thread - 1);
          partitioned_regions.emplace_back(chrom, intervals[i].start, intervals[max_i].end);
        }
      } else {
        // Too few intervals in this chromosome.
        // Create only one partition for all intervals in this chromosome.
        partitioned_regions.emplace_back(chrom, intervals[0].start, (intervals.end() - 1)->end);
      }
    }
  }
  return partitioned_regions;
}

/**
 * @brief Computes BAM and VCF features for a given region.
 * @param global_ctx Global context containing configuration and reference sequences.
 * @param worker_ctx Worker context containing BAM readers, SAM headers, and BAM indexes.
 * @param region The target region for which to compute features.
 * @param bed_regions Optional map of chromosome to vector of intervals representing BED regions.
 * @param interest_regions Optional map of chromosome to vector of intervals representing regions of interest.
 * @return A tuple containing the computed VCF features and BAM features.
 */
std::tuple<VarIdToVcfFeatures, BamRegionFeatureCollection> ComputeBamAndVcfFeaturesForRegion(
    const GlobalContext& global_ctx,
    WorkerContext& worker_ctx,
    const TargetRegion& region,
    const ChromIntervalsMap& bed_regions,
    const ChromIntervalsMap& interest_regions) {
  // This function is designed to be used by an individual worker thread as part of a larger taskflow workflow where
  // each task is getting features for a given region.

  // Compute VCF features first. Once VCF features have been computed we can use the feature set when computing BAM
  // features. This should reduce peak memory usage and lower compute time.
  // If a set of BED Regions has been passed use the intersection of the target regions and the BAM region
  const auto bam_region = Region{.chrom = region.chrom, .start = region.start, .end = region.end};
  ChromIntervalsMap target_regions;
  auto i_region = Interval{.start = region.start, .end = region.end};
  if (bed_regions.contains(region.chrom)) {
    for (const auto& interval : bed_regions.at(region.chrom)) {
      if (IntervalOverlap(interval, i_region)) {
        // Adjust interval if it falls partly outside the region
        auto interval_to_use = Interval{.start = interval.start, .end = interval.end};
        interval_to_use.start = std::max(interval_to_use.start, i_region.start);
        interval_to_use.end = std::min(interval_to_use.end, i_region.end);
        target_regions[region.chrom].emplace_back(interval_to_use);
      }
    }
  } else if (bed_regions.empty()) {
    target_regions[region.chrom].emplace_back(i_region);
  } else {
    // BED regions provided but none overlap this region; no targets
  }
  VarIdToVcfFeatures vcf_features;
  BamRegionFeatureCollection bam_features;
  const ChromPosToRefBamFeatures ref_features;
  bool use_vcf_features = true;
  if (!global_ctx.model_config.HasVcfFeatureScoringCols() || !global_ctx.model_config.use_vcf_features) {
    use_vcf_features = false;
  }
  if (!target_regions.empty()) {
    std::optional<StrUnorderedSet> skip_variants = std::nullopt;
    if (global_ctx.skip_variants_vcf.has_value()) {
      skip_variants = std::make_optional(ExtractVariantKeySet(global_ctx.skip_variants_vcf.value(), bam_region));
    }

    // Only compute VCF features for the region if there is a target entry and we have to based on workflow and scoring
    // cols
    std::optional<VarIdToVcfFeatures> passed_vcf_features = std::nullopt;
    if (use_vcf_features) {
      const VcfFeatureExtractionParams extraction_params{.is_germline_tagging = global_ctx.is_germline_tagging,
                                                         .pre_scan_variants = global_ctx.pre_scan_variants,
                                                         .pre_scan_density = global_ctx.pre_scan_density};
      const AnnotationVcfReaders annotation_readers{
          .popaf = &worker_ctx.popaf_reader, .pon = &worker_ctx.pon_reader, .catalog = &worker_ctx.catalog_reader};
      vcf_features = ExtractFeaturesForRegion(worker_ctx.vcf_feature_extraction_reader,
                                              global_ctx.genome,
                                              annotation_readers,
                                              target_regions,
                                              interest_regions,
                                              region.chrom,
                                              extraction_params);
      passed_vcf_features = vcf_features;
    }

    if (global_ctx.ref_seqs == nullptr) {
      throw error::Error("GlobalContext::ref_seqs is null — SetReferenceSequences() was not called");
    }
    bam_features =
        ComputeBamRegionFeatures(worker_ctx.alignment_readers, global_ctx.bam_feat_params, nullptr)
            .ComputeBamFeatures(bam_region, global_ctx.ref_seqs->at(region.chrom), passed_vcf_features, skip_variants);
  }
  return std::make_tuple(vcf_features, bam_features);
}

/**
 * @brief Defines a region query for a single chromosome during the parallel VCF pre-scan.
 */
struct PreScanChromTask {
  // Chromosome name
  std::string chrom;
  // htslib contig index (tid)
  s32 chrom_index;
  // 0-based start position of the region to scan
  u64 start;
  // End position (exclusive) of the region to scan
  u64 end;
};

/**
 * @brief Identifies which VCF samples to extract DP from and how they map to normal/tumor roles.
 * @details sample_indexes holds the 0-based VCF sample indexes to query. normal_pos and tumor_pos
 * are positions within sample_indexes that correspond to the normal and tumor samples respectively
 * (-1 when the role is absent). This allows ExtractDpValues to work with any number of samples
 * while still supporting the downstream ChromMedianDepth normal/tumor split.
 */
struct SampleIndexInfo {
  // 0-based VCF sample indexes to extract DP from
  vec<s32> sample_indexes;
  // Position within sample_indexes for the normal sample (-1 if absent)
  s32 normal_pos{-1};
  // Position within sample_indexes for the tumor sample (-1 if absent)
  s32 tumor_pos{-1};
};

/**
 * @brief Per-chromosome pre-scan output before merging into the global result.
 * @details Parallel workers return both partitioning entries and sample-type observations from the same
 * indexed VCF scan. The caller merges these under a mutex because both destination containers are shared.
 */
struct PreScanChromResult {
  /// Variant entries used for DP-weighted partitioning and density calculation.
  vec<VariantPreScanEntry> entries;
  /// Sample-type observations collected while scanning this chromosome region.
  SampleTypeMetricAccumulator sample_type_accumulator;
};

/**
 * @brief Mutable interval cursors for sequential VCF pre-scan.
 * @details The sequential scan advances through VCF records in coordinate order, so BED interval cursors can
 * move forward instead of rescanning intervals for every record. A single active target interval vector is
 * reused by partition filtering and sample-type detection; separate cursors are kept because those two
 * consumers mutate their cursors while applying different overlap semantics. Sample-type detection checks
 * point overlap at POS, while partition filtering checks the REF allele span. Sharing one iterator would let
 * the first check performed for a record change the candidate interval used by the second check.
 */
struct SequentialPreScanState {
  /// Active target BED intervals for the current chromosome.
  const vec<Interval>* target_intervals{};
  /// Cursor into target_intervals for partition filtering.
  vec<Interval>::const_iterator target_itr;
  /// Cursor into target_intervals for sample-type detection.
  vec<Interval>::const_iterator sample_type_target_itr;
};

vec<u32> ComputeVariantDensity(const vec<VariantPreScanEntry>& variants) {
  const auto num_variants = variants.size();
  // Each variant counts itself
  vec<u32> density(num_variants, 1);
  for (u64 i = 0; i < num_variants; ++i) {
    const auto max_pos = variants[i].pos + kPreScanDensityWindow;
    for (u64 j = i + 1; j < num_variants && variants[j].pos <= max_pos; ++j) {
      ++density[i];
      ++density[j];
    }
  }
  return density;
}

namespace {
/// Stable empty interval vector used when a chromosome has no active BED intervals.
const vec<Interval> kEmptyTargetIntervals;

/**
 * @brief Return a normalized DP value from a raw FORMAT/DP vector.
 * @details Missing, negative, and out-of-range values are normalized to zero so downstream callers can
 * use the value directly for partition weighting and sample-type tallying.
 * @param dp_values Raw FORMAT/DP values extracted from a VCF record.
 * @param sample_idx Zero-based sample index to read.
 * @return Non-negative DP value, or zero when unavailable.
 */
u32 GetDpValue(const vec<s32>& dp_values, const s32 sample_idx) {
  if (sample_idx < 0 || !std::cmp_less(sample_idx, dp_values.size())) {
    return 0;
  }
  const auto dp = dp_values.at(static_cast<size_t>(sample_idx));
  return dp > 0 ? static_cast<u32>(dp) : 0;
}

/**
 * @brief Extract DP values from raw FORMAT/DP values for each sample in the index list.
 * @param dp_values Raw FORMAT/DP values extracted from a VCF record.
 * @param sample_info Sample index information specifying which VCF samples to query.
 * @return Vector of per-sample DP values, one per entry in sample_info.sample_indexes
 */
vec<u32> ExtractDpValues(const vec<s32>& dp_values, const SampleIndexInfo& sample_info) {
  vec<u32> dps;
  dps.reserve(sample_info.sample_indexes.size());
  for (const auto idx : sample_info.sample_indexes) {
    dps.emplace_back(GetDpValue(dp_values, idx));
  }
  return dps;
}

/**
 * @brief Extract DP values from a VCF record for each sample in the index list.
 * @param record VCF record to extract DP from.
 * @param sample_info Sample index information specifying which VCF samples to query.
 * @return Vector of per-sample DP values, one per entry in sample_info.sample_indexes.
 */
vec<u32> ExtractDpValues(const io::VcfRecordPtr& record, const SampleIndexInfo& sample_info) {
  return ExtractDpValues(record->GetFormatFieldNoCheck<s32>(kFieldDp), sample_info);
}

/**
 * @brief Check whether a VCF record at the given position passes BED interval filtering.
 * Advances the BED iterator past intervals that end before the current position.
 * @param pos Zero-based VCF record position.
 * @param ref_allele Reference allele used to compute the record span.
 * @param target_intervals Sorted BED intervals for the current chromosome.
 * @param target_itr Iterator tracking the next candidate target interval; advanced in place.
 * @return true if the record passes (should be included), false if it should be skipped
 */
bool PassesBedFilter(const u64 pos,
                     const std::string_view ref_allele,
                     const vec<Interval>& target_intervals,
                     vec<Interval>::const_iterator& target_itr) {
  if (target_intervals.empty() || target_itr == target_intervals.end()) {
    return false;
  }
  while (target_itr != target_intervals.end() && target_itr->end < pos) {
    ++target_itr;
  }
  if (target_itr == target_intervals.end()) {
    return false;
  }
  return IntervalOverlap(*target_itr, pos, pos + ref_allele.length());
}

/**
 * @brief Check whether a VCF record position overlaps sample-type detection BED intervals.
 * @details Uses point-overlap semantics to match the standalone detect_sample_type app. The interval iterator
 * is advanced in place and should be reused only for sorted records on the same chromosome.
 * @param pos Zero-based VCF record position.
 * @param target_intervals Sorted BED intervals for the current chromosome.
 * @param target_itr Iterator tracking the next candidate target interval; advanced in place.
 * @return True when @p pos is inside an interval; otherwise false.
 */
bool PassesSampleTypeBedFilter(const u64 pos,
                               const vec<Interval>& target_intervals,
                               vec<Interval>::const_iterator& target_itr) {
  if (target_intervals.empty() || target_itr == target_intervals.end()) {
    return false;
  }
  while (target_itr != target_intervals.end() && target_itr->end <= pos) {
    ++target_itr;
  }
  if (target_itr == target_intervals.end()) {
    return false;
  }
  return target_itr->start <= pos && pos < target_itr->end;
}

/**
 * @brief Move completed sequential pre-scan entries into the chromosome result map.
 * @details The sequential scan accumulates one chromosome at a time. This helper centralizes the flush
 * logic used when the chromosome changes and after the final record is processed.
 * @param chrom_variants Destination map updated in place.
 * @param chrom Chromosome name for the accumulated entries.
 * @param entries Entries to move into @p chrom_variants; reset to empty after flushing.
 */
void FlushPreScanEntries(StrMap<vec<VariantPreScanEntry>>& chrom_variants,
                         const std::string& chrom,
                         vec<VariantPreScanEntry>& entries) {
  if (entries.empty()) {
    return;
  }
  chrom_variants[chrom] = std::move(entries);
  entries = {};
}

/**
 * @brief Select the target BED map used by one VCF pre-scan pass.
 * @details Sample-type detection reuses the same target BED source as filtering in filter_variants. The
 * partition BED is canonical when present; sample-type target regions are used only for callers that collect
 * sample-type metrics without partition BED filtering.
 * @param has_bed Whether partition BED filtering is enabled.
 * @param bed_regions Partition BED intervals grouped by chromosome.
 * @param sample_type_options Optional sample-type detection settings and target regions.
 * @return Non-owning pointer to the active target interval map, or nullptr when no target BED is active.
 */
const ChromIntervalsMap* SelectActiveTargetRegions(const bool has_bed,
                                                   const ChromIntervalsMap& bed_regions,
                                                   const std::optional<SampleTypePreScanOptions>& sample_type_options) {
  if (has_bed) {
    return &bed_regions;
  }
  if (sample_type_options.has_value() && sample_type_options->bed_regions.has_value()) {
    return &sample_type_options->bed_regions.value();
  }
  return nullptr;
}

/**
 * @brief Reset per-chromosome BED cursors for the sequential VCF pre-scan.
 * @details Partition filtering and sample-type filtering share the active target interval vector because
 * filter_variants uses one target BED source for both. This helper initializes independent cursors when
 * the streaming scan advances to a new chromosome.
 * @param chrom Current chromosome name.
 * @param has_bed Whether partition BED filtering is enabled.
 * @param bed_regions Partition BED intervals grouped by chromosome.
 * @param sample_type_options Optional sample-type detection settings and BED intervals.
 * @param state Mutable interval vectors and cursors reset for @p chrom.
 */
void ResetChromosomePreScanState(const std::string& chrom,
                                 const bool has_bed,
                                 const ChromIntervalsMap& bed_regions,
                                 const std::optional<SampleTypePreScanOptions>& sample_type_options,
                                 SequentialPreScanState& state) {
  const auto* const active_bed_regions = SelectActiveTargetRegions(has_bed, bed_regions, sample_type_options);

  if (active_bed_regions != nullptr && active_bed_regions->contains(chrom)) {
    state.target_intervals = &active_bed_regions->at(chrom);
  } else {
    state.target_intervals = &kEmptyTargetIntervals;
  }
  state.target_itr = state.target_intervals->begin();
  state.sample_type_target_itr = state.target_intervals->begin();
}

/**
 * @brief Add one pre-scan record to sample-type metrics using pre-extracted FORMAT/DP values.
 * @details Integrated pre-scan reads FORMAT/DP for partition entries. Reusing the raw DP vector here avoids
 * a second FORMAT/DP extraction while preserving sample-type tally semantics.
 * @param record VCF record to evaluate.
 * @param sample_type_options Sample-type detection settings.
 * @param dp_values Raw FORMAT/DP values extracted from @p record.
 * @param sample_type_accumulator Accumulator updated when the record contributes.
 */
void TallyPrecomputedDpSampleTypeRecord(const io::VcfRecordPtr& record,
                                        const SampleTypePreScanOptions& sample_type_options,
                                        const vec<s32>& dp_values,
                                        SampleTypeMetricAccumulator& sample_type_accumulator) {
  const u32 dp = GetDpValue(dp_values, sample_type_options.sample_idx);
  TallySampleTypeRecordWithDp(
      record, sample_type_options.sample_idx, dp, sample_type_options.max_af, sample_type_accumulator);
}

/**
 * @brief Append one pre-scan entry for filtering partition data.
 * @details Extracts per-sample DP values, computes total DP for region weighting, and records REF length
 * for later interval overlap calculations.
 * @param record VCF record to convert.
 * @param pos Zero-based record position.
 * @param sample_info VCF sample indexes used for DP extraction.
 * @param entries Destination entry vector updated in place.
 */
void AppendPreScanEntry(const io::VcfRecordPtr& record,
                        const u64 pos,
                        const SampleIndexInfo& sample_info,
                        vec<VariantPreScanEntry>& entries) {
  const auto ref_len = static_cast<u32>(record->Allele(0).length());
  auto dps = ExtractDpValues(record, sample_info);
  u32 total = 0;
  for (const auto dp : dps) {
    total += dp;
  }
  entries.emplace_back(pos, total, std::move(dps), ref_len);
}

/**
 * @brief Append one pre-scan entry from pre-extracted FORMAT/DP values.
 * @details Used by the integrated sample-type pre-scan path so FORMAT/DP is decoded once per record and
 * then shared by partitioning and sample-type metric tallying.
 * @param record VCF record to convert.
 * @param pos Zero-based record position.
 * @param dp_values Raw FORMAT/DP values extracted from @p record.
 * @param sample_info VCF sample indexes used for DP extraction.
 * @param entries Destination entry vector updated in place.
 */
void AppendPreScanEntryFromDpValues(const io::VcfRecordPtr& record,
                                    const u64 pos,
                                    const vec<s32>& dp_values,
                                    const SampleIndexInfo& sample_info,
                                    vec<VariantPreScanEntry>& entries) {
  const auto ref_len = static_cast<u32>(record->Allele(0).length());
  auto dps = ExtractDpValues(dp_values, sample_info);
  u32 total = 0;
  for (const auto dp : dps) {
    total += dp;
  }
  entries.emplace_back(pos, total, std::move(dps), ref_len);
}

/**
 * @brief Build per-chromosome pre-scan tasks from VCF contigs and optional BED regions.
 * @details When BED regions are provided, each task spans the BED extent for that chromosome.
 * Otherwise, each task spans the full contig length (or s64 max if the length is missing).
 * @param vcf_reader Open VCF reader used to query contig indexes and lengths
 * @param bed_regions Optional BED intervals to restrict the scan range per chromosome
 * @return Vector of PreScanChromTask, one per chromosome to scan
 */
vec<PreScanChromTask> BuildPreScanTasks(io::VcfReader& vcf_reader,
                                        const std::optional<ChromIntervalsMap>& bed_regions) {
  vec<PreScanChromTask> tasks;
  const auto hdr = vcf_reader.GetHeader();
  const auto contig_indexes = vcf_reader.GetContigIndexes();
  const auto contig_lengths = hdr->GetContigLengths();

  if (bed_regions.has_value() && !bed_regions->empty()) {
    for (const auto& [chrom, intervals] : bed_regions.value()) {
      if (contig_indexes.contains(chrom) && !intervals.empty()) {
        const auto cid = contig_indexes.at(chrom);
        tasks.emplace_back(chrom, cid, intervals.begin()->start, (intervals.end() - 1)->end);
      }
    }
  } else {
    for (const auto& [chrom, cid] : contig_indexes) {
      // Use contig length from the header when available; fall back to signed-int max so
      // htslib scans the entire contig even when ##contig length metadata is missing or zero.
      // s64 max is used because VcfReader::SetRegion forwards to htslib's hts_pos_t (int64_t).
      const u64 end = contig_lengths.contains(chrom) && contig_lengths.at(chrom) > 0
                          ? contig_lengths.at(chrom)
                          : static_cast<u64>(std::numeric_limits<s64>::max());
      tasks.emplace_back(chrom, cid, 0, end);
    }
  }
  return tasks;
}

/**
 * @brief Scan a single chromosome region for variant DP and position data using an indexed VCF query.
 * @details Opens its own VcfReader so it can be called from a parallel thread without sharing state.
 * De-duplicates by position and optionally filters by BED intervals.
 * @param vcf_file Path to the input VCF file
 * @param task Chromosome region to scan (chrom, index, start, end)
 * @param sample_info Sample index information for DP extraction
 * @param has_bed Whether BED filtering is active
 * @param bed_regions BED intervals keyed by chromosome (only consulted when has_bed is true)
 * @param sample_type_options Optional settings for collecting sample-type metrics during this scan
 * @return PreScanChromResult with sorted variant entries and optional sample-type observations for the scanned region
 */
PreScanChromResult ScanChromosomeVariants(const fs::path& vcf_file,
                                          const PreScanChromTask& task,
                                          const SampleIndexInfo& sample_info,
                                          const bool has_bed,
                                          const ChromIntervalsMap& bed_regions,
                                          const std::optional<SampleTypePreScanOptions>& sample_type_options) {
  PreScanChromResult result;
  std::unordered_set<u64> seen_positions;

  try {
    io::VcfReader vcf_reader(vcf_file);
    if (!vcf_reader.SetRegion(task.chrom_index, task.start, task.end)) {
      return result;
    }

    vec<Interval>::const_iterator target_itr;
    vec<Interval>::const_iterator sample_type_target_itr;
    const auto* const active_bed_regions = SelectActiveTargetRegions(has_bed, bed_regions, sample_type_options);
    const bool has_chrom_bed = active_bed_regions != nullptr && active_bed_regions->contains(task.chrom);
    const bool has_partition_chrom_bed = has_bed && has_chrom_bed;
    const auto& target_intervals = has_chrom_bed ? active_bed_regions->at(task.chrom) : kEmptyTargetIntervals;
    target_itr = target_intervals.begin();
    // Sample-type BED filtering is intentionally a separate cursor over the same intervals: partitioning uses
    // variant spans, while sample-type detection matches the standalone detector's point-overlap behavior.
    sample_type_target_itr = target_intervals.begin();
    const bool collect_sample_type_metrics = sample_type_options.has_value();

    while (const auto& record = vcf_reader.GetNextRegionRecord(BCF_UN_ALL)) {
      if (record->Position() < 0) {
        continue;
      }
      const auto pos = static_cast<u64>(record->Position());

      vec<s32> dp_values;
      bool has_dp_values = false;
      bool tally_sample_type_record = collect_sample_type_metrics;
      if (tally_sample_type_record && sample_type_options->bed_regions.has_value()) {
        tally_sample_type_record = PassesSampleTypeBedFilter(pos, target_intervals, sample_type_target_itr);
      }
      if (tally_sample_type_record) {
        dp_values = record->GetFormatFieldNoCheck<s32>(kFieldDp);
        has_dp_values = true;
        // Tally before de-duplicating pre-scan variant positions so integrated detection matches the
        // standalone detector, which evaluates every VCF record that overlaps the target.
        TallyPrecomputedDpSampleTypeRecord(record, *sample_type_options, dp_values, result.sample_type_accumulator);
      }

      if (has_partition_chrom_bed) {
        const bool overlaps_partition_bed = PassesBedFilter(pos, record->Allele(0), target_intervals, target_itr);
        if (!overlaps_partition_bed) {
          continue;
        }
      }

      if (!seen_positions.insert(pos).second) {
        continue;
      }

      if (has_dp_values) {
        AppendPreScanEntryFromDpValues(record, pos, dp_values, sample_info, result.entries);
      } else {
        AppendPreScanEntry(record, pos, sample_info, result.entries);
      }
    }
  } catch (const std::exception& e) {
    throw error::Error("Error in VCF pre-scan for chromosome '{}': {}", task.chrom, e.what());
  }

  std::ranges::sort(result.entries,
                    [](const VariantPreScanEntry& a, const VariantPreScanEntry& b) { return a.pos < b.pos; });
  return result;
}

/**
 * @brief Check whether a chromosome pre-scan worker produced data to merge.
 * @param result Worker-local scan result.
 * @return True when the worker produced partition entries or sample-type observations.
 */
bool HasPreScanChromResult(const PreScanChromResult& result) {
  return !result.entries.empty() || result.sample_type_accumulator.error_count > 0 ||
         !result.sample_type_accumulator.dps.empty();
}

/**
 * @brief Merge one worker-local pre-scan result into the shared pre-scan accumulators.
 * @details The worker scans with thread-local readers and containers. This helper keeps the shared-state
 * update under one lock while leaving the taskflow lambda focused on per-task work.
 * @param task Chromosome scan task associated with @p task_result.
 * @param task_result Worker-local result to merge; entries and sample-type metrics are moved from it.
 * @param result_mutex Mutex protecting @p chrom_variants and @p sample_type_accumulator.
 * @param chrom_variants Shared chromosome variant map updated with partition entries.
 * @param sample_type_accumulator Shared sample-type accumulator updated with worker observations.
 */
void MergePreScanChromResult(const PreScanChromTask& task,
                             PreScanChromResult task_result,
                             std::mutex& result_mutex,
                             StrMap<vec<VariantPreScanEntry>>& chrom_variants,
                             SampleTypeMetricAccumulator& sample_type_accumulator) {
  if (!HasPreScanChromResult(task_result)) {
    return;
  }

  const std::scoped_lock lock{result_mutex};
  if (!task_result.entries.empty()) {
    chrom_variants[task.chrom] = std::move(task_result.entries);
  }
  MergeSampleTypeAccumulators(sample_type_accumulator, std::move(task_result.sample_type_accumulator));
}

/**
 * @brief Single-threaded sequential VCF scan fallback used when no index is available or only one thread is requested.
 * @details Iterates all records in file order, grouping entries by chromosome. De-duplicates by position and
 * optionally filters by BED intervals.
 * @param vcf_file Path to the input VCF file
 * @param sample_info Sample index information for DP extraction
 * @param has_bed Whether BED filtering is active
 * @param bed_regions BED intervals keyed by chromosome (only consulted when has_bed is true)
 * @param sample_type_options Optional settings for collecting sample-type metrics during this sequential scan
 * @param sample_type_accumulator Accumulator updated with sample-type observations when options are provided
 * @return Map of chromosome name to sorted vector of VariantPreScanEntry
 */
StrMap<vec<VariantPreScanEntry>> ScanVariantsSequential(
    const fs::path& vcf_file,
    const SampleIndexInfo& sample_info,
    const bool has_bed,
    const ChromIntervalsMap& bed_regions,
    const std::optional<SampleTypePreScanOptions>& sample_type_options,
    SampleTypeMetricAccumulator& sample_type_accumulator) {
  StrMap<vec<VariantPreScanEntry>> chrom_variants;
  const io::VcfReader vcf_reader(vcf_file);
  std::string prev_chrom;
  vec<VariantPreScanEntry> entries;
  std::unordered_set<u64> seen_positions;

  SequentialPreScanState state;
  const bool collect_sample_type_metrics = sample_type_options.has_value();

  while (const auto& record = vcf_reader.GetNextRecord()) {
    if (record->Position() < 0) {
      continue;
    }
    const std::string& chrom = record->Chromosome();
    const auto pos = static_cast<u64>(record->Position());

    if (prev_chrom != chrom) {
      FlushPreScanEntries(chrom_variants, prev_chrom, entries);
      seen_positions.clear();
      prev_chrom = chrom;
      ResetChromosomePreScanState(chrom, has_bed, bed_regions, sample_type_options, state);
    }

    vec<s32> dp_values;
    bool has_dp_values = false;
    bool tally_sample_type_record = collect_sample_type_metrics;
    if (tally_sample_type_record && sample_type_options->bed_regions.has_value()) {
      tally_sample_type_record = PassesSampleTypeBedFilter(pos, *state.target_intervals, state.sample_type_target_itr);
    }
    if (tally_sample_type_record) {
      dp_values = record->GetFormatFieldNoCheck<s32>(kFieldDp);
      has_dp_values = true;
      TallyPrecomputedDpSampleTypeRecord(record, *sample_type_options, dp_values, sample_type_accumulator);
    }

    if (has_bed) {
      const bool overlaps_partition_bed =
          PassesBedFilter(pos, record->Allele(0), *state.target_intervals, state.target_itr);
      if (!overlaps_partition_bed) {
        continue;
      }
    }

    if (!seen_positions.insert(pos).second) {
      continue;
    }

    if (has_dp_values) {
      AppendPreScanEntryFromDpValues(record, pos, dp_values, sample_info, entries);
    } else {
      AppendPreScanEntry(record, pos, sample_info, entries);
    }
  }
  FlushPreScanEntries(chrom_variants, prev_chrom, entries);
  return chrom_variants;
}

/**
 * @brief Compute per-chromosome median depths from collected variant DP values.
 * @param chrom_variants Map of chromosome name to pre-scan entries containing per-variant DP values
 * @param sample_info Sample index information with normal/tumor position mappings
 * @return ChromMedianDepth with normal and tumor median DP per chromosome
 */
ChromMedianDepth ComputeChromMedianDepths(const StrMap<vec<VariantPreScanEntry>>& chrom_variants,
                                          const SampleIndexInfo& sample_info) {
  ChromMedianDepth chrom_median_dp;
  for (const auto& [chrom, entries] : chrom_variants) {
    vec<u32> normal_dps;
    vec<u32> tumor_dps;
    normal_dps.reserve(entries.size());
    tumor_dps.reserve(entries.size());
    for (const auto& entry : entries) {
      if (sample_info.normal_pos >= 0 && std::cmp_less(sample_info.normal_pos, entry.sample_dps.size())) {
        normal_dps.emplace_back(entry.sample_dps.at(static_cast<size_t>(sample_info.normal_pos)));
      }
      if (sample_info.tumor_pos >= 0 && std::cmp_less(sample_info.tumor_pos, entry.sample_dps.size())) {
        tumor_dps.emplace_back(entry.sample_dps.at(static_cast<size_t>(sample_info.tumor_pos)));
      }
    }
    if (!normal_dps.empty()) {
      chrom_median_dp.normal[chrom] = math::Median(normal_dps);
    }
    if (!tumor_dps.empty()) {
      chrom_median_dp.tumor[chrom] = math::Median(tumor_dps);
    }
  }
  return chrom_median_dp;
}

void InitializePreScanSampleInfo(const io::VcfHeaderPtr& hdr, SampleIndexInfo& sample_info) {
  if (!hdr->HasFormatField(kFieldDp)) {
    Logging::Warn("VCF file does not have FORMAT field '{}'; pre-scan DP values will be zero", kFieldDp);
  }

  const auto tn_sample_idx = GetTumorNormalSampleIndexes(hdr);
  if (tn_sample_idx.has_value()) {
    // Tumor-normal VCF: extract DP from both samples
    sample_info.sample_indexes = {tn_sample_idx->normal_sample_idx, tn_sample_idx->tumor_sample_idx};
    sample_info.normal_pos = 0;
    sample_info.tumor_pos = 1;
    return;
  }

  // Single-sample VCF: extract DP from sample 0 as normal
  sample_info.sample_indexes = {0};
  sample_info.normal_pos = 0;
}

bool RequiresSingleThreadedPreScan(io::VcfReader& vcf_reader, const fs::path& vcf_file, const size_t threads) {
  // A VCF with no records has nothing to scan, and planning per-chromosome tasks from its index would
  // fail: an index over a record-less file covers no sequence, which `GetContigIndexes` rejects.
  return threads <= 1 || !vcf_reader.HasIndex() || !VcfHasRecords(vcf_file);
}
}  // namespace

VcfPreScanResult ParallelVcfPreScan(const fs::path& vcf_file,
                                    const size_t threads,
                                    std::optional<ChromIntervalsMap> bed_regions) {
  return ParallelVcfPreScan(vcf_file,
                            {
                                .threads = threads,
                                .partition_bed_regions = std::move(bed_regions),
                                .sample_type = std::nullopt,
                            });
}

VcfPreScanResult ParallelVcfPreScan(const fs::path& vcf_file, const VcfPreScanOptions& options) {
  Logging::Info("Starting parallel VCF pre-scan of {}", vcf_file);

  // Determine sample indexes
  SampleIndexInfo sample_info;
  io::VcfReader vcf_reader(vcf_file);
  InitializePreScanSampleInfo(vcf_reader.GetHeader(), sample_info);
  const bool use_single_thread = RequiresSingleThreadedPreScan(vcf_reader, vcf_file, options.threads);
  const vec<PreScanChromTask> tasks =
      use_single_thread ? vec<PreScanChromTask>{} : BuildPreScanTasks(vcf_reader, options.partition_bed_regions);

  const bool has_bed = options.partition_bed_regions.has_value() && !options.partition_bed_regions->empty();
  const ChromIntervalsMap empty_bed_regions;
  const auto& bed_regions = has_bed ? options.partition_bed_regions.value() : empty_bed_regions;

  StrMap<vec<VariantPreScanEntry>> chrom_variants;
  SampleTypeMetricAccumulator metric_accumulator;

  if (use_single_thread) {
    chrom_variants =
        ScanVariantsSequential(vcf_file, sample_info, has_bed, bed_regions, options.sample_type, metric_accumulator);
  } else {
    Logging::Info("Pre-scanning VCF with {} threads across {} chromosomes", options.threads, tasks.size());
    std::mutex result_mutex;
    tf::Executor executor(options.threads);
    tf::Taskflow flow;
    const auto scan_task =
        [&vcf_file, &sample_info, has_bed, &bed_regions, &options, &result_mutex, &chrom_variants, &metric_accumulator](
            const PreScanChromTask& task) {
          auto task_result =
              ScanChromosomeVariants(vcf_file, task, sample_info, has_bed, bed_regions, options.sample_type);
          MergePreScanChromResult(task, std::move(task_result), result_mutex, chrom_variants, metric_accumulator);
        };
    flow.for_each(tasks.begin(), tasks.end(), scan_task).name("Parallel VCF pre-scan");
    executor.run(flow).get();
  }

  auto chrom_median_dp = ComputeChromMedianDepths(chrom_variants, sample_info);

  StrMap<vec<u32>> chrom_variant_density;
  for (const auto& [chrom, entries] : chrom_variants) {
    chrom_variant_density[chrom] = ComputeVariantDensity(entries);
  }

  Logging::Info("Pre-scan complete: {} chromosomes with variants", chrom_variants.size());

  std::optional<SampleTypeMetrics> sample_type_metrics;
  if (options.sample_type.has_value()) {
    sample_type_metrics = FinalizeSampleTypeMetrics(
        metric_accumulator, options.sample_type->total_length, vcf_file, options.sample_type->max_af);
  }

  return VcfPreScanResult{
      .chrom_variants = std::move(chrom_variants),
      .chrom_median_dp = std::move(chrom_median_dp),
      .chrom_variant_density = std::move(chrom_variant_density),
      .sample_type_metrics = sample_type_metrics,
  };
}

namespace {
/**
 * @brief Filter pre-scan entries to only those overlapping BED intervals for a single chromosome.
 * @details Uses an advancing-iterator pattern over sorted intervals and sorted entries for efficiency.
 * @param entries Pre-scan entries for one chromosome (sorted by position)
 * @param intervals Sorted BED intervals for the same chromosome
 * @return Entries that overlap at least one BED interval
 */
vec<VariantPreScanEntry> FilterEntriesByBed(const vec<VariantPreScanEntry>& entries, const vec<Interval>& intervals) {
  vec<VariantPreScanEntry> filtered;
  auto itr = intervals.begin();
  for (const auto& entry : entries) {
    const u64 var_end = entry.pos + entry.ref_len;
    while (itr != intervals.end() && itr->end < entry.pos) {
      ++itr;
    }
    if (itr != intervals.end() && IntervalOverlap(*itr, entry.pos, var_end)) {
      filtered.emplace_back(entry);
    }
  }
  return filtered;
}

/**
 * @brief Partition a single chromosome's entries into DP-weighted regions.
 * @param chrom Chromosome name
 * @param entries Pre-scan entries for this chromosome (non-empty)
 * @param threads Number of target regions
 * @param left_pad Left padding for the first region boundary
 * @param right_pad Right padding for the last region boundary
 * @param out Output vector to append regions to
 */
void PartitionChromByDp(const std::string& chrom,
                        const vec<VariantPreScanEntry>& entries,
                        const size_t threads,
                        const u32 left_pad,
                        const u32 right_pad,
                        vec<TargetRegion>& out) {
  if (entries.size() <= threads) {
    const u64 start = entries.front().pos >= left_pad ? entries.front().pos - left_pad : 0;
    const u64 end = entries.back().pos + entries.back().ref_len + right_pad;
    out.emplace_back(chrom, start, end);
    return;
  }

  u64 total_weight = 0;
  for (const auto& entry : entries) {
    total_weight += std::max(entry.total_dp, static_cast<u32>(1));
  }

  const u64 min_target = std::max(total_weight / threads, static_cast<u64>(1));
  u64 accumulated_weight = 0;
  size_t region_start_idx = 0;
  // Maximum reference end (pos + ref_len) of any entry in the current region block. A region
  // boundary is placed at the next entry's position, so cutting while the next entry starts
  // before this end would split a deletion's reference span across two regions and orphan the
  // downstream wildcard-ALT (*) records that depend on it. Track the furthest span reach so
  // the cut can be deferred until it is safe.
  u64 region_max_ref_end = 0;

  for (size_t i = 0; i < entries.size(); ++i) {
    accumulated_weight += std::max(entries[i].total_dp, static_cast<u32>(1));
    region_max_ref_end = std::max(region_max_ref_end, entries[i].pos + entries[i].ref_len);

    const bool is_last = (i == entries.size() - 1);
    // Defer an intermediate cut while the next entry falls strictly inside the current block's
    // deletion span, so deletion-spanned downstream variants stay in the same region.
    const bool span_safe = is_last || entries[i + 1].pos >= region_max_ref_end;
    if ((accumulated_weight >= min_target && span_safe) || is_last) {
      const u64 start = region_start_idx == 0 && entries[region_start_idx].pos >= left_pad
                            ? entries[region_start_idx].pos - left_pad
                            : entries[region_start_idx].pos;

      // For the last region, extend past the final variant. For intermediate regions,
      // set end = next region's first variant position to form non-overlapping half-open
      // intervals [start, end).
      const u64 end = is_last ? region_max_ref_end + right_pad : entries[i + 1].pos;

      out.emplace_back(chrom, start, end);
      region_start_idx = i + 1;
      accumulated_weight = 0;
      region_max_ref_end = 0;
    }
  }
}
}  // namespace

vec<TargetRegion> PartitionVcfRegionsByDp(const VcfPreScanResult& pre_scan,
                                          const size_t threads,
                                          const u32 left_pad,
                                          const u32 right_pad,
                                          const std::optional<ChromIntervalsMap>& bed_regions) {
  vec<TargetRegion> partitioned_regions;
  const bool has_bed = bed_regions.has_value() && !bed_regions->empty();

  for (const auto& [chrom, all_entries] : pre_scan.chrom_variants) {
    if (all_entries.empty()) {
      continue;
    }

    if (has_bed) {
      if (!bed_regions->contains(chrom)) {
        continue;
      }
      const auto filtered = FilterEntriesByBed(all_entries, bed_regions->at(chrom));
      if (!filtered.empty()) {
        PartitionChromByDp(chrom, filtered, threads, left_pad, right_pad, partitioned_regions);
      }
    } else {
      PartitionChromByDp(chrom, all_entries, threads, left_pad, right_pad, partitioned_regions);
    }
  }

  return partitioned_regions;
}

}  // namespace xoos::svc
