#include "io/output-path.h"

#include <iterator>

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/util/file-functions.h>

#include "core/read-collapser-options.h"

namespace xoos::read_collapser {

namespace {

void ValidateFileDoesNotExist(const fs::path& path, const bool overwrite) {
  if (!overwrite && file::FileExists(path)) {
    throw error::Error("Output file already exists: {}. Use --overwrite to replace.", path);
  }
}

vec<fs::path> ExpectedMarkdupOutputPaths(const ReadCollapserOptions& options) {
  vec<fs::path> paths;
  const bool mark_supplementary = options.markdup_options.mark_supplementary_alignments;
  // Per-thread BAMs are indexed unless they are merged away without a second pass over them.
  const bool write_index = mark_supplementary || !options.merge_output;
  for (size_t i = 0; i < options.threads; ++i) {
    const auto bam_path = MarkdupRegionBamPath(options, i);
    paths.emplace_back(bam_path);
    if (write_index) {
      paths.emplace_back(bam_path.string() + ".bai");
    }
    if (mark_supplementary) {
      // Phase 1 writes to temp files that Phase 2 consumes and then removes.
      const auto temp_bam_path = MarkdupRegionTempBamPath(options, i);
      paths.emplace_back(temp_bam_path);
      paths.emplace_back(temp_bam_path.string() + ".bai");
    }
  }
  if (options.merge_output) {
    const auto merged_bam_path = MarkdupMergedBamPath(options);
    paths.emplace_back(merged_bam_path);
    paths.emplace_back(merged_bam_path.string() + ".bai");
  }
  return paths;
}

vec<fs::path> ExpectedConsensusOutputPaths(const ReadCollapserOptions& options) {
  vec<fs::path> paths;
  for (size_t i = 0; i < options.threads; ++i) {
    paths.emplace_back(ConsensusFastqPath(options, i));
    if (options.consensus_options.output_cluster_bam) {
      // Cluster BAMs are not indexed because consensus may emit reads out of order.
      paths.emplace_back(ConsensusClusterBamPath(options, i));
    }
  }
  return paths;
}

}  // namespace

fs::path MarkdupRegionBamPath(const ReadCollapserOptions& options, const size_t region_index) {
  return file::PrefixedOutputPath(
      options.output_dir, options.output_prefix, fmt::format("output.{:04}.bam", region_index));
}

fs::path MarkdupRegionTempBamPath(const ReadCollapserOptions& options, const size_t region_index) {
  return file::PrefixedOutputPath(
      options.output_dir, options.output_prefix, fmt::format("output.{:04}.tmp.bam", region_index));
}

fs::path MarkdupMergedBamPath(const ReadCollapserOptions& options) {
  return file::PrefixedOutputPath(options.output_dir, options.output_prefix, "output.bam");
}

fs::path ConsensusFastqPath(const ReadCollapserOptions& options, const size_t worker_index) {
  return file::PrefixedOutputPath(
      options.output_dir, options.output_prefix, fmt::format("output.{:04}.fastq.gz", worker_index));
}

fs::path ConsensusClusterBamPath(const ReadCollapserOptions& options, const size_t worker_index) {
  return file::PrefixedOutputPath(
      options.output_dir, options.output_prefix, fmt::format("output.{:04}.bam", worker_index));
}

fs::path MetricsOutputPath(const ReadCollapserOptions& options, const std::string& default_filename) {
  // QC/statistics files live under output-dir/metrics/ per the CLI best-practices output hierarchy,
  // matching read_collapser_rescue_secondary. Deliverable BAMs/FASTQs stay in the output-dir root.
  return file::PrefixedOutputPath(options.output_dir / "metrics", options.output_prefix, default_filename);
}

vec<fs::path> ExpectedOutputPaths(const ReadCollapserOptions& options, const ReadCollapserMode mode) {
  // Metrics files are always produced, with fixed base names (optionally prefixed).
  vec<fs::path> paths{
      MetricsOutputPath(options, "summary_stats.tsv"),
      MetricsOutputPath(options, "cluster_size_distribution_summary.tsv"),
      MetricsOutputPath(options, "cluster_size_distributions.tsv"),
  };
  // `--threads` is resolved to a concrete count during parsing, so it bounds the number of per-thread files. Markdup
  // may produce fewer if the input yields fewer super regions than threads, which only makes this check conservative.
  auto mode_paths = mode == ReadCollapserMode::kMarkDuplicate ? ExpectedMarkdupOutputPaths(options)
                                                              : ExpectedConsensusOutputPaths(options);
  paths.insert(paths.end(), std::make_move_iterator(mode_paths.begin()), std::make_move_iterator(mode_paths.end()));
  return paths;
}

void ValidateOutputFilesDoNotExist(const ReadCollapserOptions& options, const ReadCollapserMode mode) {
  if (!fs::exists(options.output_dir)) {
    return;
  }
  if (!fs::is_directory(options.output_dir)) {
    throw error::Error("Output directory exists and is not a directory: {}", options.output_dir);
  }
  if (options.overwrite) {
    return;
  }
  // Only the files this run would write are checked. The output directory defaults to the current working directory,
  // so scanning every entry would both reject unrelated files and cost a full directory listing at startup.
  for (const auto& path : ExpectedOutputPaths(options, mode)) {
    ValidateFileDoesNotExist(path, false);
  }
}

}  // namespace xoos::read_collapser
