#pragma once

#include <xoos/io/metadata-util.h>

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include "bloom-filter/interleaved-bloom-filter.h"
#include "metrics/simplex-metrics.h"
#include "sequence/matcher/match-info.h"
#include "sequence/strand/strand-classifier.h"
#include "task/sink-worker.h"

namespace xoos::demux {

namespace fs = std::filesystem;

enum class ReadLengthMode {
  kFullOnly,
  kAll,
};

struct DemuxAndTrimParam {
  /// @brief A vector of input file paths containing read records to be processed.
  std::vector<fs::path> inputs;
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
  /// @brief The path to the bundle containing the possible adapter designs.
  fs::path adapter_design_bundle{"/resources/adapter-design-bundle.zip"};
  /// @brief The adapter design name from the bundle to be used.
  // TODO: Consider making the downstream of this non-optional (It is nice to see the default value in the help)
  // i.e. LoadAdapterDesign is expecting optional. Perhaps address this during a redesign of adapter design bundle
  std::string adapter_design_name{"SBX-D"};
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
  /// @brief The maximum error rate allowed in intra-molecular consensus (pre-trim) duplex data
  float max_error_rate_percent{10.0f};
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

  // Control parameters not directly exposed to users, useful because this struct remains const throughout runtime
  /// @brief If true, the data will be processed as methylation data
  bool methylation{false};

  /// @brief Comments to be added to output files
  io::Comments comments{};
};

namespace detail {

using SidPool = std::unordered_map<uint, Barcode>;

SidPool LoadSampleSheet(const BarcodePool& sid_pool, const fs::path& samplesheet);

BarcodePool LoadSampleSheet(const fs::path& samplesheet);

SidPool LoadSidPool(const BarcodePool& sid_pool);

}  // namespace detail

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
