<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Demux will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- Added dual-index simplex adapter support and the example `I2X-S-DX` design.
  - Added `--dual-sid-mode`, which allows partial assignment when a SID is unique in a pair.
- Added new simplex, duplex, and run-level metrics; see the Changed section for a summary.
- Added a quality score for duplex bases at topologically ambiguous indels, controlled by `--ambiguous-topology-quality-mode`.
  - The new quality is counted in the `ambiguous_duplex_bases` metric.
- Added `--quality-score-overrides` to set the four duplex quality characters.

### Changed

- CLI options have changed for cross-module consistency.
  - `demux` changes:
    - `--input` → `--reads`
    - `--input-files-list` → `--reads-file-list`
    - `--out-dir` → `--output-dir`
    - `--writing-threads-per-sample` → `--writers-per-sample`
    - `--worker-threads-per-input` → `--workers-per-reads-file`
    - `--min-trimmed-read-len` → `--min-trimmed-read-length`
    - `--output-failed-reads` → `--write-failed-reads`
    - `--length-distribution-report-max` → `--max-read-length-bin`
    - `--min-read-len` → `--min-read-length`
    - `--stop-after-min-concordant-duplex-bases` → `--min-concordant-duplex-bases-for-early-stop`
    - `--strand-detect` → `--strand-index`
    - `--max-error-rate-percent` → `--max-pairwise-error-rate`; use decimal values from `0` to `1`.
    - The `-b`, `-p`, and `-n` short aliases are removed. Some legacy long-option aliases remain available for temporary compatibility but are deprecated.
  - `demux_strand_index` changes:
    - `--input` → `--reference`
    - `--output` → `--output-strand-index`
    - All short aliases are removed.
- Additional CLI-related changes:
  - Removed `--suppress-simplex-qual-override`. Simplex input qualities are preserved by default; use `--override-simplex-quality` to apply fixed quality values.
  - The default for `--discordant-sid-mode` is now `discard-all`, which discards all discordant SID reads. Use `discard-tied` to retain the previous clear-winner behavior.
  - The `--adapter-design-name` option is now required; it has no default.
  - Added `--quality-score-overrides` for runtime configuration of base-call quality characters.
  - Added `--ambiguous-topology-quality-mode`; the default `smear` mode emits the fourth quality score for ambiguous duplex bases.
- Output changes:
  - Read-name bitflag now populates the threaded-state bits (bit 7 `startsThreaded`, bit 6 `endsThreaded`) from the upstream `molecular_traces_read_state` RDB annotation, instead of always zeroing the upper bits.
  - Samples that receive no output reads now produce a per-sample directory containing a single empty FASTQ, providing a valid output for every sample-sheet entry.
  - `run_stats.tsv` now uses `value` instead of `count` and adds a `denominator` column.
  - Replaced single histogram TSVs with run-level and per-sample distributions under `metrics/read_length_distr/`.
  - `sample_stats.tsv` now includes percentage rows.
  - Simplex metrics have been overhauled to better match duplex metrics and other requirements:
    - New run-level and sample-level metrics:
      - `total_bases`
      - `trimmed_bases`
      - `mean_passing_read_length`
      - `sid_concordant_reads`
      - `filtered_bases`
      - `full_only_discarded_reads`
      - `umi_failed_reads`
    - New run-level metrics:
      - `raw_bases`
      - `unassigned_bases`
      - `failed_assigned_bases`
      - `trimmed_bases`
      - `passing_reads`
      - `too_long_reads`
    - Renamed `full_reads` → `full_simplex_reads`
    - Renamed `partial_reads` → `partial_simplex_reads`
    - Removed `index_hopping_reads`
    - Simplex `unassigned_reads`, `too_short_trimmed_reads`, and `failed_assigned_reads` reflect classification after SID assignment matching duplex behavior.
  - Read-length histogram output was overhauled:
    - Output to subdirectory `metrics/read_length_distr/`.
    - File name prefixes are now: `raw-input`, `raw-failed`, `no-hairpin`, `passing-read`, `end-adapter-position`
- Adapter changes:
  - Standardized default HPSIDs to a canonical orientation for every adapter design that uses HPSIDs. Use canonical SID orientation in sample sheets. Reverse-complemented SIDs are accepted for legacy compatibility.
  - YSU:
    - SID and UMI detection is now independent; SID-only matches are trimmed to the SID edge.
    - At least one end must contain both a SID and UMI; otherwise, the read is discarded as `umi_failed_reads`.
  - Switched duplex pairwise alignment from infix to prefix mode, improving alignment speed and increasing duplex read yield by better handling repeated sequences within a read.
- Restructured GitBook to multi-page layout.

### Fixed

- Simplex fixes:
  - `too-short-after-trim` is now determined after SID assignment, matching duplex.
  - Reads that match no SID are now classified as `unassigned_reads` only.
  - Fixed an issue where `--read-length-mode all-split` put incorrect SIDs in the read FASTQ names of the partial fraction.
  - Simplex SID searches no longer use `--min-trimmed-read-length` as an early search limit. The option remains a post-trim filter, preventing retained adapter sequence and classifying short assigned reads as `too_short_trimmed_reads` instead of `unassigned_reads`.

### Removed

- Removed the redundant reverse-complemented HPSID FASTA files for YS, YS-NEW, and YSU from the adapter bundle.

### Known issues

- **Sample-count limit:** More than 256 samples cannot be processed in one run. Split larger sheets into batches and rerun on `raw_failed`; assignments may favor earlier batches. Support for more than 256 samples per run is planned.
- **Existing YS behavior identified during 1.2.0 validation:** In a no-adapter negative control with 105 SIDs, YS wrote 0.49% of reads to sample output, compared with 0.03% for the slower YS-NEW design.
  More SIDs increase opportunities for a non-specific single-SID match; smaller pools reduce but do not eliminate this behavior.
  Consider YS-NEW if the YS rate is too high. For YS-NEW, increase `--min-score` and revalidate assignment sensitivity; use `--read-length-mode full-only` to discard partial reads when they are not needed.
- **Partial simplex reads falsely classified as short full-length reads:** In long-read libraries using YS-NEW, dual-index simplex, or other experimental generic non-UMI simplex designs, partial reads can be overtrimmed after a low-confidence paired SID match triggers a non-specific adapter match.
  Review `full_passing_reads` and `passing_reads` in `metrics/read_length_distr/total_trimmed_simplex_passing_read_length_distr.tsv`.
  For custom affected designs whose SID entries use `max_edit_distance` `2` (for example, YS-NEW), lowering that value to `1` can reduce non-specific matches but reduces sensitivity to genuine SID errors; `0` may be too stringent.
  Alternatively, filter full reads by the FASTQ read bitflag and length, or use `--read-length-mode all-split` and filter only the `full` output.
  These mitigations can reduce the impact but cannot completely prevent the issue; a code fix is required.
  Current YS uses a separate non-probabilistic scoring and trimming methodology and remains available as a fallback.
- **Long-standing SBX-D long-read trimming limitation:** Present since the first public Demux release and identified during 1.2.0 validation, a long SBX-D consensus can be trimmed at an incorrect boundary when the double-read heuristic selects an unrelated start-adapter match.
  SBX-D read-length distributions can show an artificial inflection around 430 bases.
  Review long-read trim boundaries before downstream analysis; no configuration workaround is available in 1.2.0.
- **Custom adapter support added in 1.2.0:** Custom designs have limitations with fixed terminal sequences shorter than 8 bases, nonstandard duplex linker or stem structures, and exceptionally complex simplex fixed regions.
  Use bundled designs where possible and validate custom designs before analysis; bundled adapter designs are not affected.

## [1.1.0]

### Changed

- Changed output FASTQ filename format to `<sample_name>-<part_index>.fastq[.gz|.zst]` to align with the Everest naming convention. Part indices are 1-based, zero-padded (6 digits by default), contiguous, and assigned per sample in input file order.
- Input files are now sorted lexicographically before processing to ensure deterministic part index assignment.
- Migrated TSV/histogram metadata from legacy `#`-prefixed comments to structured `##RocheCommandLine=<ID=...,Version="...",CommandLine="...">` header format.
- Using `NA` instead of `0` for metrics not applicable to duplex UMI adapter types and for strand metrics when strand detection is disabled.
- Sample names are determined from prefixes instead of directory names, supporting nested output directories.
- CLI `--help` output now groups adapter-specific options under "Simplex Options" and "Duplex Options" headings for easier discovery.
- Changed discordant SID behavior of simplex reads to no longer partial trim or trim both SIDs when discordant.

### Added

- Added `--read-length-mode all-split` option that writes full and partial simplex reads into separate subdirectories (`<sample>/full/` and `<sample>/partial/`).
- Added `--suppress-simplex-qual-override` flag. By default, base quality scores for simplex adapter reads (Ys, YsuTe, YsuSl, Simplex, etc.) are overridden with a fixed simplex quality value. When this flag is set, original quality scores are preserved.
- Added new simplex adapter type (`YS-NEW`) for improved SID assignment accuracy over the existing `YS` adapter.
- Added new `SIMPLEX-10X` adapter type for 10x-linked simplex libraries. Only the runway and SID are trimmed. No trimming on the 3'-end is applied.
- Added `--min-score` CLI parameter to set a minimum log-odds score threshold for adapter matching (default `30`, experimental).
- Added metrics raw_bases, unassigned_bases, failed_assigned_bases, and trimmed_bases to duplex metrics to help with potential troubleshooting.
- Added sample sheet validation to reject the reserved sample name `raw_failed`.
- Added `file_provenance.tsv` in each sample output directory, mapping each output FASTQ file to its source input file.
- Added optional `prefix` and `suffix` fields to LUT definitions in the adapter design manifest. Affixes are prepended/appended to FASTA barcode cores during LUT construction, keeping FASTA files as pure barcode sequences while encoding adapter linker context in the manifest.
- Added TSV support for `--sample-sheet`. The delimiter (comma or tab) is auto-detected.
- Added `--worker-threads-per-input` CLI option to control how many input files are processed concurrently (replaces the hardcoded `kCoresPerDataset = 4` heuristic). Default `4` preserves existing behavior.
- Added optional `bait_5p` and `bait_3p` fields in adapter design manifest to specify expected 5' and 3' bait sequences for each adapter type.
- Added `--discordant-sid-mode` CLI option to control handling of reads with discordant 5'/3' SIDs. Three modes: `discard-tied` (discards reads where both sides are equally likely), `discard-all` (discards all discordant reads), `keep` (keeps all, assigns to winner). In all modes, only the winning side's trim boundary is applied and 5' wins ties.
- Added `sid_discordant_discarded_reads` and `failed_assigned_reads` metrics to simplex run metrics.

### Fixed

- Fixed `IsPartialRead` for YSU adapter: reads with one SID and zero UMIs are now correctly classified as partial. Previously they were neither full nor partial, breaking the `assigned = full + partial` invariant.
- Improved duplex hairpin detection: fixed out-of-bounds Bitap search window near read ends, corrected SID search anchoring for indel-containing loop matches, and scaled search windows to loop length. Recovers ~0.04% additional assigned reads with zero throughput regression.
- Fixed non-deterministic `file_sid` assignment caused by `unordered_map` iteration in sample sheet parsing. Sample sheet entries are now processed in row order, making `file_sid` values stable across runs.
- Fixed `raw_failed` FASTQ files being written to the output root instead of the `raw_failed/` subdirectory.
- Fixed simplex metrics counting: `preassignment_passing_reads` now includes reads filtered after trimming, and `too_short_reads` vs `too_short_trimmed_reads` are tracked separately.
- Fixed YSU false positive UMI matches by extending the LUT match region from 12-mer to 15-mer (8-mer core + `CG` prefix + `GAGTT` suffix) and gating SID trimming on UMI confirmation. Eliminates spurious SID-only trims on reads without a valid UMI.
- Fixed `SeqMatcher::FindNextBarcode` failing to find exact barcode matches when shorter approximate matches existed at adjacent start positions. The function now tries all candidate lengths at each start position and selects the best by log-likelihood score (exact matches win immediately). For 3' scanning, an extra look-back of `max_len - min_len` positions ensures longer exact matches are not shadowed by shorter approximate hits at higher start positions.

### Removed

- Removed unused `stem-5p.fa` and `stem-3p.fa` from `YSU` adapter design bundle.

## [1.0.0]

### Added

- Adapter trimming and sample demultiplexing for SBX reads.
- Adapter designs: SBX-D (duplex), SBX-DM (TAPS+ methylation with XM tag output), SBX-FAST (simplex), YS.
- SID identification and UMI extraction from adapter sequences.
- YC tag annotation (v1.5) with base-type quality scores.
- Duplex and simplex read processing with configurable read-length filtering.
- FASTQ input validation for sequence/quality length consistency and valid characters.
- Compressed output via gzip (default) or zstd.
- Per-sample output directories with configurable writing threads.
- Run-level and sample-level metrics output (`run_stats.tsv`, `sample_stats.tsv`).
- Strand index generation via `demux_strand_index`.

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.0.0...1.1.0
[1.0.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.0.0
