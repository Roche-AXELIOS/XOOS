#pragma once

#include <cstddef>
#include <string>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

#include "core/read-collapser-options.h"

namespace xoos::read_collapser {

// Output filename helpers. These are the single source of truth for how output files are named, so
// that the writer sites and ExpectedOutputPaths (used for pre-run overwrite checks) always agree.
// Each joins options.output_prefix to the default filename with an underscore when the prefix is
// non-empty (e.g. prefix "SAMPLE" -> "SAMPLE_output.bam"), matching xoos::file::PrefixedOutputPath.

// Per-super-region markdup BAM, e.g. "<prefix_>output.0003.bam".
fs::path MarkdupRegionBamPath(const ReadCollapserOptions& options, size_t region_index);

// Per-super-region markdup temp BAM used when --mark-supplementary-alignments is enabled.
fs::path MarkdupRegionTempBamPath(const ReadCollapserOptions& options, size_t region_index);

// Merged markdup BAM produced when --merge-output is set, e.g. "<prefix_>output.bam".
fs::path MarkdupMergedBamPath(const ReadCollapserOptions& options);

// Per-worker consensus FASTQ, e.g. "<prefix_>output.0003.fastq.gz".
fs::path ConsensusFastqPath(const ReadCollapserOptions& options, size_t worker_index);

// Per-worker consensus cluster BAM, e.g. "<prefix_>output.0003.bam".
fs::path ConsensusClusterBamPath(const ReadCollapserOptions& options, size_t worker_index);

// Metrics TSV path for @p default_filename (e.g. "summary_stats.tsv").
fs::path MetricsOutputPath(const ReadCollapserOptions& options, const std::string& default_filename);

// Returns the exact set of output paths that a run with @p options in @p mode may produce.
vec<fs::path> ExpectedOutputPaths(const ReadCollapserOptions& options, ReadCollapserMode mode);

// Checks that none of the output files this run would produce already exist. Only the paths returned by
// `ExpectedOutputPaths` are checked, so unrelated files in the output directory are left alone. This matters because
// the output directory defaults to the current working directory, which is routinely non-empty.
// Throws xoos::error::Error if the output path is not a directory, or if output files exist and overwrite is false.
void ValidateOutputFilesDoNotExist(const ReadCollapserOptions& options, ReadCollapserMode mode);

}  // namespace xoos::read_collapser
