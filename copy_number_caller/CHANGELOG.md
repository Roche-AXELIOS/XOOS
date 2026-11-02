<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Copy Number Caller will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- Added the `tumor-only-wgs` subcommand for somatic copy-number calling when no matched normal sample is available.
- Added automatic selection of bundled seed segments, mappability tracks, and region blocklists for WGS workflows. Explicit resource options continue to override these defaults, and `germline-wgs` can detect or accept a library type to select the appropriate interval size and mappability track.
- Added `--output-prefix` to prepend a sample prefix to output filenames and `--overwrite` to allow replacement of existing outputs.
- Added a low-purity fallback to whole-genome-doubling (WGD) QC. When predicted purity is below 0.35, the diploid-bias solution is rejected in favour of a uniform copy-number prior if the assessable autosomal genome fraction (by length) with loss of heterozygosity reaches the `--min-loh-fraction-for-wgd` threshold (default 0.35).

### Changed

- Changed WGS `--sex` handling to support `auto-detect`, `male`, `female`, and `skip-prediction` strategies (`unknown` remains a deprecated alias), and report resolved values as `male`, `female`, or `unknown` instead of `M`, `F`, or `N`.
- **Breaking:** Reorganized default outputs into `metrics/`, `tracks/`, and `supplemental/` directories while keeping primary callsets in the output root.
- Standardized output conventions: missing values use `NA`, and VCF/TSV metadata uses the XOOS `##` + PascalCase keys.

### Fixed

- Fixed hg38 seed segments so they exclude the true centromere.
- Fixed a segmentation bug where a tiny BAF (DH) segment contained within a single LogR observation could trigger a runaway merge that swallowed downstream segments, dropping a legitimate boundary such as a CN-LOH to balanced transition. Boundaries are now re-resolved after each merge, so these transitions are preserved. Affects the `germline-wgs` and `tumor-normal-wgs` subcommands.

### Removed

- Removed the bundled hg19 seed-segment and mappability resources.

### Known issues

- Performance may be limited in tumors with low SCNV signal, such as those with low tumor purity and/or low SCNV burden. See [Low SCNV signal tumors](docs/troubleshooting.md#low-scnv-signal-tumors) in the Troubleshooting guide for analysis guidance.
- Purity and ploidy prediction remains challenging in tumors with genome doubling because of model identifiability: multiple purity/ploidy combinations may explain the data equally well. See the [Whole genome doubled tumors](docs/troubleshooting.md#whole-genome-doubled-tumors) in the Troubleshooting guide for analysis guidance.
- Germline CNV calling from Giraffe pangenome-aligned BAMs may miss common duplications unless a rescue step is applied. See the [Call germline copy numbers](docs/usage-germline-wgs.md) page for details.

## [1.1.0]

### Added

- Genomic copy number variation detection by analyzing coverage patterns in DNA sequencing alignments.
- Germline WGS workflow (`germline-wgs`) for detecting copy number variants in normal samples.
- Somatic tumor-normal WGS workflow (`tumor-normal-wgs`) for characterizing allele-specific somatic CNVs in matched tumor-normal pairs.
- Segmentation workflow (`segmentation`) for standalone genome segmentation.
- Purity/ploidy estimation with grid search and whole genome doubling (WGD) QC.
- Somatic CNV prediction via `predict-somatic-cnv` subcommand.
- Interval coverage calculation with GC-content normalization.
- Intermediate output files in standard BED format with metadata headers.
- VCF and SEG output with version and command line metadata for traceability.
- Mirrored BAF filter and low interval density filter for purity/ploidy estimation.
- Decoy interval filtering via `hg38_decoy_intervals.bed` to reduce false positive duplications.
- Low-coverage early exit: when median coverage is below 10x, the module writes valid header-only output files.
- Graceful exit on VCF data-quality issues in `tumor-normal-wgs` (low het-SNP fraction, zero variants after filtering).
- Warning and skip for calls in pseudo-autosomal regions (PAR) if chrY PAR alignments are detected.
- Configurable blocklist BED for excluding problematic regions.

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.1.0
