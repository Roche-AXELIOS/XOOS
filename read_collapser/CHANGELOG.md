<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Read Collapser will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- Added `--output-prefix` to `markdup` and `consensus`, and a `##Sample=<...>` header to metrics TSVs for sample-aware outputs.

### Changed

- Replaced `--include-per-base-read-support-tags` and `--include-per-base-majority-count-tags` with `--include-consensus-aux-tags`. The new option emits standard typed consensus auxiliary tags, with disagreeing-base error counts replacing majority counts.
- Renamed `read_collapser_rescue_secondary --prefix` to `--output-prefix` and changed the prefix separator from `.` to `_`; the old flag is no longer accepted.
- Moved markdup and consensus metrics into `output-dir/metrics/` and clarified empty-result metrics by distinguishing measured `0` from unavailable `NA` values.
- Renamed CLI parameters for cross-module consistency. The old names are no longer accepted:
  - `--bam-input` → `--bam`
  - `--bed-input` → `--bed`
  - `--output-cluster-bam` → `--write-cluster-bam`
  - `--consensus-threshold` → `--min-consensus-fraction`
  - `--consensus-gap-threshold` → `--min-consensus-gap-fraction`
  - `--min-mapq` → `--min-mapping-quality`
  - `--skip-softclips` → `--skip-softclip-consensus`
- The default output directory for `markdup` and `consensus` is now the current directory (`.`), matching `read_collapser_rescue_secondary`.
- `read_collapser_rescue_secondary` will now reject BAM files that have already been processed by `read_collapser_rescue_secondary`, preventing accidental re-processing of the same data as secondary rescue is not an idempotent operation.

### Fixed

- Fixed duplex reads with missing UMIs being treated as full-length, so partial reads can be assigned to the correct UMI cluster.
- Fixed a bug where `markdup` would add a PG line to the SAM header with a non-unique ID if a PG line with the same ID already existed. A unique suffix is now appended to the new ID to avoid collision.

## [1.1.0]

### Added

- Added `--overwrite` flag to `markdup` and `consensus` subcommands. Both now fail fast if output files already exist unless `--overwrite` is specified.
- Added `wgs-simplex` preset to consensus.
- Added `--ignore-read-name-parsing-errors` as a hidden CLI option for processing BAMs with non-standard read names.
- Added `--mark-supplementary-alignments` option to `markdup`. When enabled, supplementary alignments whose primary alignment was marked as duplicate are also flagged. A new metric `duplicate_supplementary_alignments` is reported.

### Changed

- Soft-clip consensus is now enabled by default. The `--include-softclips` flag has been replaced with `--skip-softclips` to opt out.
- Duplicate marking now prefers full-length reads over partial reads when selecting the representative non-duplicate read in a cluster.
- The `wgs-simplex` markdup preset has been updated: removed `--cluster-by-strand` and `--cluster-by-umi`, changed `--wiggle-room` from `0` to `2`, changed `--wiggle-room-partial` from `2` to `0`.
- Read names are now parsed for all reads to determine full/partial read status via adapter bitfield, not only in the UMI clustering path.
- Metric TSV files now use a single `##RocheCommandLine=<...>` metadata line instead of the previous two-line format.
- Downsampling seed changed from a fixed value to a hash of the cluster ID to avoid bias. This may cause minor differences in consensus sequences for downsampled clusters.
- Renamed metrics in `summary_stats.tsv` to use consistent terminology distinguishing records, alignments, and reads (e.g. `input_reads` → `input_records`, `unclustered_supplementary_reads` → `unclustered_supplementary_alignments`).
- Unmapped reads missing UMIs are now properly discarded and reflected in the `discarded_missing_umi_records` metric.

### Fixed

- Fixed merging of optional per-base tag vectors when concatenating soft-clip consensus results.
- Fixed off-by-one in `TrimEnds` that discarded a single surviving base when `trimmed_start` equaled `trimmed_end`.
- Fixed a crash when generating soft-clip consensus for large clusters where not all reads overlap.

## [1.0.0]

### Added

- Duplicate marking via `markdup` subcommand: position-based and UMI-aware clustering of aligned reads, with configurable region padding and BED file input.
- Consensus generation via `consensus` subcommand: error-corrected consensus reads from duplicate clusters, with per-base depth and majority count tags.
- Presets for common workflows: `wgs-duplex`, `wgs-duplex-mrd`, `wgs-duplex-cfdna`, `wgs-simplex`.
- Configurable minimum cluster size, strand-specific cluster size thresholds, and simplex depth filtering.
- Partial read clustering via `--cluster-partials`.
- Merged output mode via `--merge-output` for `markdup`.
- Summary statistics output with detailed read-level metrics (clustering input, unclustered, discarded, and consensus counts).

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.0.0...1.1.0
[1.0.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.0.0
