<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Alignment Metrics will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- `read_metrics_summary.tsv` now reports `unique_aligned_bases` and `unmapped_soft_clipped_bases` (and their `_passing_filter` counterparts), derived from each read's primary and supplementary (SA tag) alignments.
- Added `--stratify-by-label` to stratify metrics by region class in a single run. When set, the 4th column of `--bed-input` is read as a comma-separated list of labels (up to 16), overlapping regions are split at label boundaries, and metrics are attributed independently to each label and written to per-label subdirectories of the output directory (`OUTPUT_DIR/<label>/{accuracy,coverage,read_metrics}`). Without the flag, the output layout is unchanged.
- Added the `merge_and_label_bed` tool, which builds the labeled BED consumed by `--stratify-by-label` from one or more plain BED files paired with labels.
- `--bed-input` and `merge_and_label_bed --bed` now accept gzip- or bgzip-compressed BED files named with a `.gz` extension, in addition to uncompressed BED files.
- Added `--overwrite` flag to error out when output files already exist, unless explicitly overridden.
- Added `--output-prefix` parameter to prepend a prefix to all output filenames.
- The tool now logs a warning for zero-read inputs and produces valid empty metrics output.
- Added the `gcbias` subcommand to compute GC bias metrics independently of coverage, accuracy, and read metrics.
It requires a BAM file and a reference FASTA, and writes *gc_bias/gc_bias.tsv* and *gc_bias/gc_bias_summary.tsv* to the output directory.
- Added `--enable-gc-bias-metrics` to the `coverage` and `all` subcommands to include GC bias output alongside other metrics from those subcommands.
- Added `alignment_records` to `read_metrics_summary.tsv` to report primary, secondary, and supplementary alignment records.
- The tool now warns when `--min-base-quality` or `--min-base-type` is incompatible with the detected duplex, simplex, or intermolecular-consensus dataset type.
- Added `read_length_histograms.tsv` and `read_length_summary.tsv` reporting read length including soft-clipped bases. The read length excluding soft-clipped bases metrics are now written to `read_length_excluding_soft_clips_histograms.tsv` and `read_length_excluding_soft_clips_summary.tsv`.
- Added `mapped_*` read length columns (`mapped_read_count`, `mapped_partial_read_count`, `mapped_full_read_count` and the corresponding summary columns) that count all mapped primary alignments before read-level filters are applied. These appear alongside the existing `post_filter_*` columns, which continue to count the subset of reads passing the mapping quality and flag filters.

### Changed

- Reorganized CLI help menu into standard sections.
- Renamed CLI parameters for consistency across the suite:
  - `--bam-input` → `--bam`
  - `--bed-input` → `--bed`
  - `--min-mapq` → `--min-mapping-quality`
  - `--min-baseq` → `--min-base-quality`
  - `--trim-leading-bases` → `--trim-start-length`
  - `--trim-trailing-bases` → `--trim-end-length`
  - `--disable-base-type-decoding` → `--skip-base-type-decoding`
  - `--exclude-uncovered-positions` → `--skip-uncovered-positions`
  - `--enable-te-metrics` → `--enable-target-enrichment-metrics`
  - `--hp-*`/`--enable-hp-*` → `--homopolymer-*`/`--enable-homopolymer-*`
- Removed all single-letter shorthands except `-h`, `-v`, and `-l`.
- Moved `--enable-target-enrichment-metrics` into its own `Target Enrichment Metrics Options` help section, shared across the `coverage`, `read`, and `all` subcommands, and restored it for `read` (previously only available for `coverage` and `all`).
- Adopted shared deep input-file validators for BAM, FASTA, and BED inputs.
- `aligned_bases` now includes aligned bases from supplementary alignments (from the SA tag) that overlap a region of interest, and may double-count, so its percentage is reported as `NA`; `unique_aligned_bases` is reported as a percentage of `total_bases` instead.
- Restructured GitBook to multi-page layout.
- All floating point metric values now use two decimal places, except error percentages, which use six decimal places. Percentage values no longer include a `%` suffix when the metric's row or column name already identifies the value as a percentage.
- Coverage and accuracy metrics now exclude reference positions with non-ACGT bases from calculations.
- Added `--include-ambiguous-reference-bases` option to `coverage`, `accuracy`, and `all` subcommands. When enabled, only N/n reference positions are excluded; IUPAC ambiguous bases (R, Y, S, W, K, M, B, D, H, V) are included as valid positions and query bases are not counted as errors when consistent with the ambiguity code.
- Reference input `--reference` is now required for `coverage` subcommand.
- Removed the duplex-only `concordant_duplex_coverage` metric from coverage outputs. All dataset types now report `any_coverage` and `post_filter_coverage`.

### Fixed

- Fixed zero-denominator Phred scores in `base_level_accuracy_summary.tsv` producing Q93 instead of `NA` for zero-read inputs.

### Known issues

- The allele fraction filter `--max-alt-allele-fraction` is intended to exclude likely variants from being counted as errors in the metrics.
This filter tends to overestimate the accuracy of long homopolymers because it may cause systematic errors to be excluded. This is particularly noticeable for simplex data where long homopolymers may paradoxically appear to be more accurate than shorter ones. To mitigate this, users should take into consideration the potential impact of the allele fraction filter when interpreting homopolymer accuracy metrics and disable the allele fraction filter by setting `--max-alt-allele-fraction 1.0` if necessary. This issue does not affect cases where the allele fraction filter is not applied. Whenever possible, use a BED file that excludes known variant sites (e.g. the high confidence regions for GIAB samples) and disable the allele fraction filter for the most reliable homopolymer accuracy assessment.

## [1.1.0]

### Changed

- Replaced legacy metadata in TSV outputs with structured `##RocheCommandLine=<...>` metadata format. All TSV files now include tool name, version, and full command line in a single metadata line.
- Zero coverage positions are now included in the calculation of mean, median, percentiles, and ratio metrics when a BED file is provided. This produces more accurate statistics when metrics are restricted to specific regions.

## [1.0.0]

### Added

- Read-level metrics: total reads, aligned reads, unmapped reads, forward/reverse strand counts, duplicate counts, and base counts (total, aligned, soft-clipped, unmapped).
- Coverage metrics: per-base and binned coverage depth, GC-bias statistics.
- Accuracy metrics: mismatch, insertion, and deletion rates with homopolymer-aware breakdowns.
- Homopolymer accuracy metrics with per-length error rates and spanning/effective read counts.
- Target enrichment metrics via `--enable-te-metrics`: on-target rate, mean/median coverage, fold-80 base penalty, uniformity.
- Subcommands: `read`, `coverage`, `accuracy`, `all`.
- YC tag base-type decoding for SBX duplex and simplex chemistries, with `--disable-base-type-decoding` for standard BAM files.
- Configurable read filtering via `--exclude-flags`, `--min-mapq`, and trimming of leading/trailing bases.

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.0.0...1.1.0
[1.0.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.0.0
