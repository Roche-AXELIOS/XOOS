#pragma once

#include <condition_variable>
#include <memory>
#include <mutex>
#include <optional>
#include <queue>
#include <string>
#include <vector>

#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

#include "compute-bam-features/alignment-reader.h"
#include "compute-bam-features/compute-bam-region-features.h"
#include "core/config.h"
#include "core/score-calculator.h"
#include "file_type_inference/sample-type-metrics.h"
#include "util/region-util.h"

namespace xoos::svc {

/**
 * Synopsis:
 * This file contains utility functions and data structures for parallel feature computation for filtering variants.
 */

// Window size (upstream and downstream) for computing variant density during the pre-scan pass.
// A 201-bp window means 100 bp upstream + the variant position + 100 bp downstream.
static constexpr u64 kPreScanDensityWindow{100};

/**
 * @brief Per-variant metadata collected during the VCF pre-scan pass.
 * @details Stores position, depth, and allele length for each unique-position variant. Used for
 * DP-weighted region partitioning and pre-computed variant density.
 */
struct VariantPreScanEntry {
  // 0-based position
  u64 pos{};
  // sum of DP across all samples (used for region weighting)
  u32 total_dp{};
  // per-sample DP values, indexed parallel to the sample index list
  vec<u32> sample_dps;
  // length of REF allele (for interval computation)
  u32 ref_len{};
};

/**
 * @brief Result of the unified parallel VCF pre-scan.
 * @details Replaces the two separate sequential VCF passes (GetChromosomeMedianDepth + PartitionVcfRegions)
 * with a single parallel pass that collects DP values, variant positions, and variant density.
 */
struct VcfPreScanResult {
  // Per-chromosome variant data, sorted by position
  StrMap<vec<VariantPreScanEntry>> chrom_variants;
  // Computed chromosome median depths (from the collected DP values)
  ChromMedianDepth chrom_median_dp;
  // Pre-computed variant density per chromosome, indexed parallel to chrom_variants
  StrMap<vec<u32>> chrom_variant_density;
  // Optional sample-type metrics collected during pre-scan
  std::optional<SampleTypeMetrics> sample_type_metrics;
};

/**
 * @brief Optional sample-type detection inputs for VCF pre-scan.
 * @details When present in VcfPreScanOptions, the pre-scan collects SampleTypeMetrics while it is already
 * reading VCF records for filtering. The BED regions are independent from partitioning regions because
 * sample-type detection uses point-overlap semantics and its own denominator length.
 */
struct SampleTypePreScanOptions {
  /// Zero-based VCF sample index evaluated for sample-type metrics.
  s32 sample_idx{};
  /// Maximum allele frequency for observations counted as low-frequency errors.
  f64 max_af{kDefaultMaxAfForSampleTypeDetection};
  /// Denominator length used to finalize the sample-type phred metric.
  u64 total_length{kDefaultWholeGenomeLengthForSampleTypeDetection};
  /// Optional target intervals restricting records evaluated for sample-type detection.
  std::optional<ChromIntervalsMap> bed_regions;
};

/**
 * @brief Configuration for the unified VCF pre-scan.
 * @details Bundles existing partitioning options with optional sample-type detection so callers can choose
 * one VCF pass that serves both filtering setup and deferred sample-type inference.
 */
struct VcfPreScanOptions {
  /// Requested number of pre-scan threads.
  size_t threads{1};
  /// Optional BED regions used for variant partitioning and filtering setup.
  std::optional<ChromIntervalsMap> partition_bed_regions;
  /// Optional sample-type detection settings; absent when pre-scan should not collect metrics.
  std::optional<SampleTypePreScanOptions> sample_type;
};

/**
 * @brief Optional external-VCF input paths used to annotate features.
 * @details Bundles the POPAF, Panel of Normals (PON), and catalog VCF paths into a single parameter to
 * avoid passing three adjacent same-type arguments to the WorkerContext constructor.
 */
struct AnnotationVcfPaths {
  std::optional<std::string> popaf;
  std::optional<std::string> pon;
  std::optional<std::string> catalog;
};

/**
 * @brief Manages per-thread mutable state for variant processing.
 *
 * Each thread owns its own WorkerContext to avoid concurrent access. Contains VCF readers,
 * score calculators, and alignment readers that are reused across regions processed by the
 * same thread.
 */
struct WorkerContext {
  io::VcfReader filter_variants_reader;
  io::VcfReader vcf_feature_extraction_reader;
  std::optional<io::VcfReader> popaf_reader;
  std::optional<io::VcfReader> pon_reader;
  std::optional<io::VcfReader> catalog_reader;
  vec<ScoreCalculator> calculators;
  vec<AlignmentReader> alignment_readers;

  /**
   * @brief Construct a WorkerContext with readers and calculators for the given inputs.
   * @param vcf_file Path to the input VCF file.
   * @param annotation_vcfs Optional POPAF/PON/catalog VCF paths used to build the annotation readers.
   * @param bam_inputs Paths to the input BAM files.
   * @param alignment_reader_cache Shared cache for BAM index data across threads.
   */
  WorkerContext(const std::string& vcf_file,
                const AnnotationVcfPaths& annotation_vcfs,
                const std::vector<fs::path>& bam_inputs,
                AlignmentReaderCache& alignment_reader_cache);
};

/**
 * @brief Per-region output data populated by a filter task. Contains VCF records and SHAP value rows
 * for a single region.
 */
struct RegionResult {
  vec<io::VcfRecordPtr> out_records;
  vec<vec<std::string>> shap_value_rows;
  vec<vec<std::string>> snv_shap_value_rows;
  vec<vec<std::string>> indel_shap_value_rows;
};

/**
 * @brief Item pushed to the writer queue by filter tasks. Carries the region's output data and its
 * index so the writer can associate it with the correct per-region temp file.
 */
struct RegionOutput {
  size_t region_index{};
  RegionResult result;
};

/**
 * @brief Thread-safe queue for passing completed region results from filter tasks to the writer thread.
 *
 * Supports a shutdown sentinel (empty optional) to signal the writer that no more items will arrive.
 * After Shutdown() is called, Push() silently drops items and Pop() returns nullopt once the queue
 * drains, allowing both producer and consumer threads to terminate cleanly.
 */
class RegionOutputQueue {
 public:
  /**
   * @brief Enqueue a completed region result.
   *
   * After Shutdown() has been called, items are silently dropped to prevent unbounded
   * memory growth when the writer has already failed or finished.
   *
   * @param item Region output containing the region index and its filtered records/SHAP rows.
   */
  void Push(RegionOutput item) {
    {
      const std::lock_guard lock(_mutex);
      if (_shutdown) {
        return;
      }
      _queue.emplace(std::move(item));
    }
    _cv.notify_one();
  }

  /**
   * @brief Signal that no more items should be accepted or waited for.
   *
   * Wakes any thread blocked in Pop(). Subsequent Push() calls are no-ops.
   * Pop() will return nullopt once all previously enqueued items have been consumed.
   */
  void Shutdown() {
    {
      const std::lock_guard lock(_mutex);
      _shutdown = true;
    }
    _cv.notify_one();
  }

  /**
   * @brief Block until an item is available, then dequeue and return it.
   *
   * Returns nullopt when the queue has been shut down and all remaining items have been consumed.
   *
   * @return The next RegionOutput, or nullopt if the queue is shut down and empty.
   */
  std::optional<RegionOutput> Pop() {
    std::unique_lock lock(_mutex);
    _cv.wait(lock, [this]() { return !_queue.empty() || _shutdown; });
    if (_queue.empty()) {
      return std::nullopt;
    }
    auto item = std::move(_queue.front());
    _queue.pop();
    return item;
  }

  /**
   * @brief Check whether the queue has been shut down.
   * @return True if Shutdown() has been called.
   */
  bool IsShutdown() const {
    const std::lock_guard lock(_mutex);
    return _shutdown;
  }

 private:
  mutable std::mutex _mutex;
  std::condition_variable _cv;
  std::queue<RegionOutput> _queue;
  bool _shutdown = false;
};

// Forward declaration — defined in util/locked-tsv-writer.h.
class LockedTsvWriter;

/**
 * @brief Groups the optional SHAP value TSV writers passed to the writer thread.
 *
 * Each writer is nullable — when absent, the corresponding SHAP output is skipped.
 * Writers are shared pointers because they may be referenced by both the writer thread
 * and the main thread for header/footer operations.
 */
struct ShapWriters {
  /// Combined SNV + InDel SHAP values.
  std::shared_ptr<LockedTsvWriter> combined;
  /// SNV-only SHAP values.
  std::shared_ptr<LockedTsvWriter> snv;
  /// InDel-only SHAP values.
  std::shared_ptr<LockedTsvWriter> indel;
};

/**
 * @brief The global state required to filter variants, this state is immutable and shared across all threads.
 */
struct GlobalContext {
  fs::path genome;
  ComputeBamFeaturesParams bam_feat_params{};
  const StrMap<std::string>* ref_seqs{nullptr};
  ChromIntervalsMap bed_regions;
  ChromIntervalsMap interest_regions;
  // Regions whose start positions receive the germline-mitochondria `blocklist` FILTER tag.
  ChromIntervalsMap blocklist_regions;
  SVCConfig model_config;
  io::VcfHeaderPtr hdr;
  vec<io::InfoFieldMetadata> vcf_info_metadata;
  vec<io::FormatFieldMetadata> vcf_fmt_metadata;
  std::optional<s32> vcf_normal_index;
  std::optional<s32> vcf_tumor_index;
  ChromMedianDepth normalize_targets{};
  bool normalize_scoring_features{};
  sex_predict::Sex sex{};
  std::string chr_x_name;
  std::string chr_y_name;
  vec<Interval> chr_x_par;
  vec<Interval> chr_y_par;
  bool phased{};
  std::optional<StrMap<vec<Interval>>> force_calls;
  std::optional<StrUnorderedSet> hotspots;
  std::optional<StrUnorderedSet> block_list;
  f32 min_af{};
  // Minimum depth threshold for the germline-mitochondria workflow's `low_dp` filter.
  u32 min_dp{};
  // AF threshold at or above which a germline-mitochondria variant is genotyped as homozygous ALT (1/1).
  f32 min_hom_alt_af{};
  f32 weighted_counts_threshold{};
  f32 hotspot_weighted_counts_threshold{};
  f32 ml_threshold{};
  f32 snv_min_ml_score{};
  f32 indel_min_ml_score{};
  f32 hotspot_ml_threshold{};
  f32 germline_fail_ml_threshold{};
  f32 min_phased_allele_freq{};
  f32 max_phased_allele_freq{};
  u32 min_alt_counts{};
  std::optional<fs::path> skip_variants_vcf;
  u32 min_tumor_support{};
  // tumor-only-wgs: minimum single-sample alt support hard-filter threshold.
  u32 min_support{};
  u32 max_normal_support{};
  f32 min_tumor_af{};
  // Relaxed thresholds applied to variants flagged `in_catalog` (--catalog-vcf).
  // When unset, the corresponding standard threshold is used (no relaxation).
  std::optional<f32> catalog_min_tumor_af;
  std::optional<f32> catalog_snv_min_ml_score;
  std::optional<f32> catalog_indel_min_ml_score;
  f32 min_dp_ratio{};
  u32 max_indel_size{};
  bool is_germline_tagging{};
  // Pre-computed variant density from the VCF pre-scan, keyed by chromosome.
  // When set, ExtractFeaturesForRegion uses these values instead of computing density per-region.
  const StrMap<vec<VariantPreScanEntry>>* pre_scan_variants{nullptr};
  const StrMap<vec<u32>>* pre_scan_density{nullptr};
};

/**
 * @brief Compute BAM and VCF features for all variants in a single region.
 *
 * Reads variants from the VCF, extracts alignment-level features from the BAM, and
 * computes VCF-level features. Results are returned as a tuple of per-variant VCF features
 * and per-region BAM features.
 *
 * @param global_ctx Immutable shared state (reference genome, model config, header, etc.).
 * @param worker_ctx Per-thread mutable state (VCF readers, alignment readers, score calculators).
 * @param region Genomic region to process.
 * @param bed_regions BED intervals used for depth normalization.
 * @param interest_regions BED intervals defining regions of interest for filtering.
 * @return Tuple of (per-variant VCF features, per-region BAM feature collection).
 */
std::tuple<VarIdToVcfFeatures, BamRegionFeatureCollection> ComputeBamAndVcfFeaturesForRegion(
    const GlobalContext& global_ctx,
    WorkerContext& worker_ctx,
    const TargetRegion& region,
    const ChromIntervalsMap& bed_regions,
    const ChromIntervalsMap& interest_regions);

/**
 * @brief Partition a VCF file into roughly equal regions for parallel processing.
 *
 * Scans the VCF to count variants per contig and splits them into regions sized for
 * the given thread count. Optionally restricts to BED intervals.
 *
 * @param vcf_file Path to the input VCF file.
 * @param threads Number of threads (determines target region count).
 * @param bed_regions Optional BED intervals to restrict partitioning to.
 * @return Vector of target regions covering all variants in the VCF.
 */
vec<TargetRegion> PartitionVcfRegions(const fs::path& vcf_file,
                                      size_t threads,
                                      std::optional<ChromIntervalsMap> bed_regions);

/**
 * @brief Perform a single parallel pass through a VCF file to collect per-variant DP values, variant positions,
 * and variant density.
 * @details Replaces the two separate sequential VCF passes (GetChromosomeMedianDepth + ExtractVariantIntervals)
 * with a single parallel pass. Each chromosome is processed in its own thread using indexed region queries.
 * Falls back to single-threaded if the VCF has no index or only one thread is requested.
 * @param vcf_file Path to the input VCF file
 * @param threads Number of threads for parallel processing
 * @param bed_regions Optional BED regions to restrict which variants are included
 * @return VcfPreScanResult containing per-chromosome variant data, median DPs, and variant density
 */
VcfPreScanResult ParallelVcfPreScan(const fs::path& vcf_file,
                                    size_t threads,
                                    std::optional<ChromIntervalsMap> bed_regions);

/**
 * @brief Perform a configurable VCF pre-scan in one sequential or parallel pass.
 * @details Extends the legacy pre-scan with optional sample-type metric collection. The same VCF read
 * collects per-variant DP, density inputs, partitioning inputs, and, when requested, sample-type phred
 * observations. This avoids a separate VCF pass before filtering.
 * @param vcf_file Path to the input VCF file. Indexed VCFs are scanned per contig when multiple threads are requested.
 * @param options Threading, partition BED, and optional sample-type detection settings.
 * @return VcfPreScanResult containing partition data, chromosome median depths, density values, and optional
 * sample-type metrics.
 */
VcfPreScanResult ParallelVcfPreScan(const fs::path& vcf_file, const VcfPreScanOptions& options);

/**
 * @brief Partition pre-scanned variants into regions weighted by DP for balanced parallel filtering.
 * @details Instead of splitting so each region has an equal number of variants, splits so each region has
 * approximately equal total DP weight. This ensures high-depth regions are split into smaller regions with
 * fewer variants, preventing them from becoming bottlenecks during parallel filtering.
 * @param pre_scan Pre-scan result containing per-chromosome variant data
 * @param threads Number of threads (determines target number of regions per chromosome)
 * @param left_pad Left padding for region boundaries
 * @param right_pad Right padding for region boundaries
 * @param bed_regions Optional BED regions to restrict which pre-scan variants are used for partitioning.
 *        When provided, only variants overlapping a BED interval are included in the partition.
 *        This is separate from the pre-scan BED filter so that median DP and density can be computed
 *        from the full VCF while partitioning is restricted to BED regions.
 * @return Vector of TargetRegion objects for parallel filtering
 */
vec<TargetRegion> PartitionVcfRegionsByDp(const VcfPreScanResult& pre_scan,
                                          size_t threads,
                                          u32 left_pad,
                                          u32 right_pad,
                                          const std::optional<ChromIntervalsMap>& bed_regions);

// Overload without BED filtering.
inline vec<TargetRegion> PartitionVcfRegionsByDp(const VcfPreScanResult& pre_scan,
                                                 const size_t threads,
                                                 const u32 left_pad,
                                                 const u32 right_pad) {
  return PartitionVcfRegionsByDp(pre_scan, threads, left_pad, right_pad, std::nullopt);
}

// Overload with default padding (1 bp on each side) and no BED filtering.
inline vec<TargetRegion> PartitionVcfRegionsByDp(const VcfPreScanResult& pre_scan, const size_t threads) {
  return PartitionVcfRegionsByDp(pre_scan, threads, 1, 1);
}

// Overload with default padding (1 bp on each side) and BED filtering.
inline vec<TargetRegion> PartitionVcfRegionsByDp(const VcfPreScanResult& pre_scan,
                                                 const size_t threads,
                                                 const std::optional<ChromIntervalsMap>& bed_regions) {
  return PartitionVcfRegionsByDp(pre_scan, threads, 1, 1, bed_regions);
}

/**
 * @brief Compute variant density for a sorted vector of variant positions using a sliding window.
 * @details Uses a 201-bp window (100 bp upstream + variant position + 100 bp downstream). Each variant's
 * density is initialized to 1 (itself) and incremented for each other variant within the window.
 * Operating on the full chromosome avoids boundary artifacts that occur when density is computed per-region.
 * @param variants Sorted vector of pre-scan entries for a single chromosome
 * @return Vector of density values, indexed parallel to the input variants
 */
vec<u32> ComputeVariantDensity(const vec<VariantPreScanEntry>& variants);

}  // namespace xoos::svc
