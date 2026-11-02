<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to STR Caller will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- Added `RareRepeatFXN` support for FXN-style repeats, including off-target read recruitment for BWA inputs.
- Added `--pileup-visualization-loci` to limit pileup SVGs to selected catalog loci.
- Added experimental `--call-str-region-variants` to call SNVs and indels in and around STR loci and annotate the STR VCF with `LocusSequenceVariant`.
- Added `--overwrite` to replace existing VCF/JSON outputs; without it, existing outputs fail early, and `--output-dir` must be a directory.

### Changed

- The packaged hg38 catalog is now used by default.
- `--aligner` replaces `--giraffe-input-bam` and accepts `unknown`, `bwa`, or `giraffe`; `unknown` is the default and auto-detects the aligner from the BAM header, falling back to `bwa` when no supported aligner is found.
- Renamed CLI options for consistency:
  - `--bam-input` → `--bam`
  - `--exclude-sex-output` → `--skip-sex-output`
  - `--weighted-purity-score-cutoff` → `--min-weighted-purity-score`
  - `--sp-fraction-hom-cutoff` → `--min-spanning-reads-hom-fraction`
  - `--sp-fraction-low-confidence-cutoff` → `--max-spanning-reads-low-confidence-fraction`
  - `--metrics` → `--output-metrics`
  - `--write-to-bam-unmapped` → `--write-reads-as-unmapped`
- Metrics output is now written to `<output-dir>/metrics/metrics.tsv` with the standardized tabular layout.
- `--sex` now accepts `auto-detect`, `male`, `female`, or `skip-prediction`; `--sex unknown` is no longer valid.

### Fixed

- Runs with simplex or unconfirmed-duplex input BAMs now warn and write empty outputs.
- JSON output now writes unavailable confidence intervals as `null` instead of `"NaN"`.

### Known Issues

- FMR1 coverage is insufficient for reliable STR genotyping in AXELIOS 1 sequencer software version 1.0.0 and earlier due to a known issue in the intra-consensus generation step that impacts repeat estimates.
- STR Caller supports duplex (SBX-D and SBX-FAST) input. For simplex (SBX-S) or
  unconfirmed-duplex BAMs, the caller logs a warning and writes empty output
  files.
- Non-determinism in the BWA alignment of in-repeat reads—where reads consisting entirely of a repeat motif are randomly distributed across genomic regions with long identical repeats, depending on input order or sorting—can skew repeat count estimates for loci like FMR1 when expansions exceed the read length.

## [1.1.0]

### Added

- Genotyping and detection of repeat expansions in short tandem repeats (STR) from SBX sequencing data.
- VCF and JSON output describing alternate repeat structures.
- Metrics TSV output with per-locus statistics.
- STR catalog support for defining target loci.
- Support for SBX duplex reads with lossless encoding to recover R1 and R2 reads for additional support.
- Configurable region padding, read filtering, and confidence thresholds.
- FXN locus support with configurable spanning-read fraction cut-offs for HOM and low-confidence calls.
- Support for compressed reference files.

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.2.0
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.1.0
