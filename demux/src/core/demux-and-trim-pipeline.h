#pragma once

#include <filesystem>  // NOLINT(build/c++17)
#include <string>

#include <xoos/types/vec.h>

#include "bloom-filter/interleaved-bloom-filter.h"
#include "core/dual-sid-mode.h"
#include "core/trimming-debug-mode.h"
#include "metrics/simplex-metrics.h"
#include "sequence/strand/strand-classifier.h"
#include "task/sink-worker.h"
#include "utility/alignment-util.h"

namespace xoos::demux {

namespace fs = std::filesystem;

enum class ReadLengthMode {
  kFullOnly,
  kAll,
  kAllSplit,
};

/// @brief Controls assignment of the ambiguous duplex quality score;
enum class AmbiguousTopologyQualityMode {
  kSmear,
  kNone,
  kSmearOverride,
};

struct DemuxAndTrimParam {
  /// @brief A vector of input file paths containing read records to be processed.
  vec<fs::path> inputs;
  /// @brief If provided then only these samples are demultiplexed and output files named accordingly
  std::optional<fs::path> sample_sheet;
  /// @brief The output directory where trimmed read records and metrics will be stored.
  fs::path out_dir{"."};
  /// @brief Allow overwrite of files in output rather than throwing an error if they already exist.
  bool overwrite{false};
  /// @brief The number of threads to be used for parallel processing.
  size_t threads{1};
  /// @brief The batch size of read records to be processed together.
  size_t batch_size{500};
  /// @brief Read filtering type for demultiplexing simplex reads
  ReadLengthMode read_length_mode{ReadLengthMode::kAll};
  /// @brief Controls whether dual-SID simplex reads may be assigned from only one detected SID side.
  DualSidMode dual_sid_mode{DualSidMode::kNone};
  /// @brief How to handle reads where 5' and 3' SIDs disagree.
  DiscordantSidMode discordant_sid_mode{DiscordantSidMode::kDiscardAll};
  /// @brief Assignment mode for the ambiguous duplex quality score on topologically ambiguous duplex bases.
  /// Defaults to kSmear so the 4th (ambiguous) duplex quality score is emitted by default
  AmbiguousTopologyQualityMode ambiguous_topology_quality_mode{AmbiguousTopologyQualityMode::kSmear};
  /// @brief Output quality scores (phred+33) for simplex, concordant, discordant, and ambiguous bases.
  BaseQualityScores base_quality_scores{};
  /// @brief The path to the bundle containing the possible adapter designs.
  fs::path adapter_design_bundle{"/resources/adapter-design-bundle.zip"};
  /** @brief Required adapter design name from the bundle's DSL config. */
  std::string adapter_design_name;
  /// @brief The compression level for output files, 1 = fastest, 9/19 = slowest (zlib/zstd).
  std::optional<size_t> compression_level;
  /// @brief Minimum read length of the sequence for processing to continue.
  size_t min_read_len{50};
  /// @brief Minimum trimmed read length of the sequence for processing to continue.
  size_t min_trimmed_read_len{50};
  /// @brief The maximum length to report in the read length distribution metrics.
  u32 length_distribution_report_max{1000};
  /// @brief The compression type used for the output FASTQ files.
  SinkCompressionType compression_type{SinkCompressionType::kGzip};
  /// @brief The maximum number of writing threads allocated per SID
  size_t writing_threads_per_sample{1};
  /**
   * @brief Divisor used to determine how many input files are processed concurrently.
   * Concurrent inputs = threads / worker_threads_per_input (minimum 1). Clamped to threads if larger.
   */
  size_t worker_threads_per_input{4};
  /// @brief The maximum error rate allowed in intra-molecular consensus (pre-trim) duplex data.
  f32 max_pairwise_error_rate{0.1F};
  /// @brief Minimum concordant duplex bases to process before stopping. If set to std::nullopt should not stop early.
  std::optional<u64> min_concord_dp_bases;
  /// @brief For strand classification critical value as the upper bound for false strand classification
  // 70 is 10^(-7) FPR
  double strand_critical{70};
  /// @brief Loaded from filename for interleaved Bloom Filter for strand detection mode.
  std::optional<strand::InterleavedBloomFilter> strand_detector;
  /// @brief The strand classifier for strand classification mode.
  std::optional<strand::StrandClassifier> strand_classifier;
  /// @brief Output of non-demuxed raw reads, usually for debugging purposes
  bool output_failed_reads{false};
  /// @brief When true, override base quality scores for simplex adapter reads with the configured simplex score
  bool override_simplex_qual_score{false};
  /// @brief Minimum score needed for a valid match log odds 30 is 1/(2^30) = 1/(1024^3) is roughly equal to 1/(10^9)
  // that is we have 1 in a billion chance of a false positive per site (assuming our scoring scheme reflects reality)
  // This a match is counts as 2 so this is a minimum of 15 exact base matches
  s32 min_score{30};

  // Control parameters not directly exposed to users, useful because this struct remains const throughout runtime
  /// @brief If true, the data will be processed as methylation data
  bool methylation{false};
  /// @brief Trimming debug modes: emit extra tags in output FASTQ comments for trimming evaluation.
  vec<TrimmingDebugMode> trimming_debug;

  /// @brief Tool name, version, and command line for output file metadata
  io::CommandLineInfo command_line_info{};
};

/**
 * @brief Check metrics validation and handle existing metrics files.
 *
 * Creates the metrics directory if it doesn't exist. If overwrite is not enabled
 * and any metrics files already exist, an error will be thrown.
 * @param param The demux and trim parameters
 * @throws error::Error if metrics directory is not a directory or if metrics files exist without overwrite flag
 */
void HandleMetricsFiles(const DemuxAndTrimParam& param);

/**
 * @brief Executes the demultiplexing and trimming pipeline for input reads.
 *
 * This function orchestrates a pipeline that performs demultiplexing and trimming of input reads.
 * The pipeline includes reading input read records, demultiplexing them based on the provided strategy,
 * and trimming adapters. The trimmed records are then written to an output directory along with metrics.
 */
void DemuxAndTrimPipeline(const DemuxAndTrimParam& param);

}  // namespace xoos::demux
