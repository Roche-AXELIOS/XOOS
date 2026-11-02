#pragma once

#include <optional>
#include <string>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/yc-decode/yc-decoder.h>

namespace xoos::alignment_metrics {

constexpr s32 kDefaultMaxClusterSizeBin = 50;
constexpr s32 kDefaultMaxReadLengthBin = 1'000;
constexpr s32 kDefaultMaxReferenceBufferLength = 1'000;
constexpr s32 kDefaultRegionSize = 10'000;
constexpr u32 kDefaultMinDepth = 10;
constexpr f64 kDefaultMaxAltAlleleFraction = 0.1;
constexpr u32 kDefaultMaxCoverageBin = 10'000;
constexpr u32 kDefaultMinBaseq = 30;
constexpr u8 kDefaultBaseQualityThresholdForHpMasking = 5;
constexpr auto kDefaultMinBaseType = yc_decode::BaseType::kSimplex;
constexpr auto kMaxTrimLeadingBases = 50U;
constexpr auto kMaxTrimTrailingBases = 50U;
constexpr u32 kDefaultMinMapq = 4;
constexpr u32 kMinBaseQuality = 0U;
constexpr u32 kMaxBaseQuality = 255U;
const std::vector<u64> kDefaultSummaryStatsPercentiles{10, 25, 50, 75, 90};
const std::vector<u64> kDefaultCoverageCutoffs{10, 20, 30};
constexpr f64 kDefaultHpSubsamplingFraction = 1.0;
constexpr u8 kHpMinLength = 2;
constexpr u8 kHpMaxLength = 30;
constexpr u32 kRequiredAnchorOverlap = 1;
constexpr u16 kDefaultExcludeFlag = BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FUNMAP | BAM_FDUP;
constexpr u32 kDefaultGcWindowSize = 100;
constexpr u32 kMaxGcWindowSize = 1'000;

constexpr u32 kMinClusterSize = 1U;
constexpr u32 kMaxClusterSize = 255U;

struct MetricsTypes {
  bool has_coverage_metrics{true};
  bool has_accuracy_metrics{true};
  bool has_read_metrics{true};
  bool has_gc_bias_metrics{false};
};

// Parameters for coverage metrics.
struct AlignmentMetricsOptions {
  fs::path bam_input;
  fs::path bam_index;
  fs::path reference_input;
  fs::path bed_input;
  // When true, interpret the BED 4th column as a comma-separated list of region labels and stratify
  // metrics output into one subdirectory per label. Requires bed_input.
  bool stratify_by_label{};

  fs::path out_dir;
  // Prefix prepended to all output filenames (joined with an underscore). When empty, output
  // filenames are unchanged.
  std::string output_prefix;
  // When true, allow WriteMetrics to overwrite existing output files. Without this flag, WriteMetrics
  // throws if any output file it would write already exists.
  bool overwrite{};
  vec<u64> coverage_cutoff;
  vec<u64> summary_stats_percentiles;
  bool ignore_family{};
  bool exclude_empty_positions{};
  bool disable_hp_quality_modification{};
  u8 base_quality_threshold_for_hp_masking{};

  u32 reference_padding{};
  u32 max_coverage_bin{};
  u32 min_depth{};
  u32 max_cluster_size_bin{};
  u32 max_read_length_bin{};
  u32 region_size{};
  u32 gc_scan_window_size{};
  f64 max_alt_allele_fraction{};
  u16 exclude_flags{};
  u16 min_mapq{};
  u16 min_baseq{};
  yc_decode::BaseType min_base_type{};
  bool disable_base_type_decoding{};

  size_t threads{};

  bool calculate_hp_metrics{};
  u16 min_hp_length{};
  u16 max_hp_length{};
  bool hp_allow_heterogeneous_insertions{};

  f64 hp_subsampling_fraction{kDefaultHpSubsamplingFraction};
  std::optional<int> hp_subsampling_seed{std::nullopt};

  // Base trimming options
  u16 trim_leading_bases{0};
  u16 trim_trailing_bases{0};

  bool include_ambiguous_reference_bases{false};

  bool enable_te_metrics{};

  MetricsTypes metric_types{};
  io::CommandLineInfo command_line_info{};
};

/**
 * Checks whether the given AlignmentMetricsOptions indicates that
 * homopolymer masking of read qualities is needed.
 */
bool NeedHpMasking(const AlignmentMetricsOptions& options);

/**
 * Checks whether the given AlignmentMetricsOptions indicates that
 * decoding of base types is needed.
 */
bool NeedBaseTypeDecoding(const AlignmentMetricsOptions& options);

}  // namespace xoos::alignment_metrics
