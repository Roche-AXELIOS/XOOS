<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Tumor Fraction Estimator will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- Added contamination filtering by YC-tag base type with `--min-base-type` (default `simplex`), plus `--output-prefix` and `--overwrite` options for output control.
- Added `--tumor-sample-name` (required) and `--normal-sample-name` (optional) for named-sample selection in `tumor-fraction`, with support for multi-sample, reordered, and multi-allelic VCFs.
- Added dedicated low-quality ALT filter statuses and metrics so ALT-only filtering is reported separately from read and pileup filtering.

### Changed

- Renamed CLI options for clarity and consistency:
  - `--min-baseq` → `--min-base-quality`
  - `--min-mean-baseq` → `--min-mean-base-quality`
  - `--min-mapq` → `--min-mapping-quality`
  - `--vcf` → `--query-vcf`
  - `--max-germline-allele-depth` → `--max-normal-ad`
  - `--min-somatic-vaf` → `--min-tumor-af`
- Renamed variant-level filter statuses (in the `filter` column of the counts TSV and the `filtered_variants_*` metric rows):
  - `GermlineAD` → `HighNormalAD`
  - `LowSomaticVAF` → `LowTumorAF`
  - `ObservedVAF` → `HighObservedAltVAF`
  - `OtherObservedVAF` → `HighObservedNonQueryAltVAF`
  - `AltCoverage` → `HighAltCoverage`
  - `NoiseProbeOverlappingSomatic` → `NoiseProbeOverlappingQueryVariant`
- Renamed read/base observation filter statuses (in the `filtered_reads_*` metric rows and `alt_*` counts columns):
  - `MAPQ` → `MappingQual`
  - `Mismatch` → `MismatchCount`
  - `Softclip` → `SoftclipLength`
  - `Indel` → `IndelProximity`
- Renamed variant type `Somatic` to `QueryVariant` to avoid confusion when running in contamination mode with population SNPs.
- Renamed results metrics and counts columns to use `non_query` instead of `other`:
  - `total_other_alts` → `total_non_query_alts`
  - `mean_vaf_other_alts` → `mean_vaf_non_query_alts`
  - `noise_total_other_alts` → `noise_total_non_query_alts`
  - `noise_mean_vaf_other_alts` → `noise_mean_vaf_non_query_alts`
  - counts column `other_alts_count` → `non_query_alt_count`
- Updated the default values for the following CLI options for the `tumor-fraction` subcommand:
  - `--min-coverage` default changed from 10 → 30
  - `--max-softclip-length` default changed from 0 → 50
  - `--min-tumor-af` default changed from 0.03 → 0.1
  - `--max-mismatch` default changed from 5 → 3
- Restructured output filenames for consistency:
  - `tumor_fraction.results.tsv` → `tumor_fraction_estimation_results.tsv`
  - `contamination.results.tsv` → `contamination_estimation_results.tsv`
  - `tumor_fraction.counts.tsv` / `contamination.counts.tsv` → `per_site_counts.tsv`
- Replaced the bundled included-regions BED with a GIAB high-confidence region set that excludes difficult-to-map regions and the Y chromosome, improving estimation accuracy.
- Contamination fractions below 0.0001 are now reported as `<0.0001`.

### Fixed

- Corrected included- and excluded-region matching to evaluate each variant's one-base interval, and corrected BAM `N` (reference-skip) CIGAR handling so skipped reference positions are treated as uncovered rather than deletions.

## [1.1.0]

### Added

- Tumor fraction estimation from SBX duplex sequencing data via the `tumor-fraction` subcommand.
- Contamination estimation via the `contamination` subcommand using population SNP VCFs.
- TSV metrics output with per-subcommand prefixed filenames (`tumor_fraction.results.tsv`, `contamination.results.tsv`).
- Input validation with file extension and readability checks, and fail-fast errors when output files already exist.
- Read filtering controls grouped into read-level, base-level, pileup, and noise-estimation options.
- `--save-read-features` developer option for emitting per-read feature data.

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.1.0
