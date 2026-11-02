<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Small Variant Caller will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- Add full `tumor-only-wgs` support across SVC submodules, including bundled pre-trained BWA-aligned FFPE and cell-line models.
- Add SBX Simplex (SBX-S) support to germline WGS feature computation, model training, and variant filtering through the `germline-wgs-simplex` profile and bundled BWA and Giraffe SNV and indel models.
- Add the `germline-mitochondria` workflow to call mitochondrial small variants from GATK Mutect2 output using allele-frequency and depth thresholds instead of an ML model. The workflow adds subcommands to `filter_variants`, `compute_bam_features`, and `compute_vcf_features`; `filter_variants` adds `--min-af`, `--min-dp`, `--min-hom-alt-af`, and `--blocklist-bed`, the `low_af` and `low_dp` FILTER tags, and drops spanning-deletion (`*`) alleles.
- Add tumor-normal sample-type detection through `--sample-type auto-detect` and the `detect_sample_type` utility, with optional detection controls.
- Add VCF-based QC metrics output to `filter_variants`. Every workflow now writes a `metrics/metrics.tsv` summary (input/output/PASS variant counts, PASS genotype breakdown, het/hom and Ti/Tv ratios, and a per-FILTER breakdown), with a `qc_status` field that reports `fail` when there are no output or no PASS variants. Adds the `--output-metrics` CLI option (default `metrics/metrics.tsv`, relative to `--output-dir`).
- Add catalog-aware annotation and filtering with the `in_catalog` feature, `--catalog-vcf`, and catalog-specific tumor AF and ML score thresholds for `tumor-normal-wgs`.
- Bundle a pre-built catalog of recurrent TCGA somatic hotspots (GRCh38) at `resources/tcga_hotspots.vcf.gz`.
- Add duplex-derived BAM features and use the AF-based genotype call to rescue otherwise failed germline diploid reconciliation, including multi-allelic calls.
- Add `phase_variants`, a utility for phasing somatic variant calls after filtering with `filter_variants tumor-normal-wgs` and `filter_variants tumor-only-wgs`. Variants passing filtering are assigned phase groups, with an option to combine complex variants and MNVs into one call.
- Add `features_to_vcf`, a utility that converts BAM feature TSV rows with `chrom`, `pos`, `ref`, and `alt` columns into an indexed single-sample VCF with synthetic heterozygous genotype fields for downstream SVC feature workflows.

### Changed

- Consolidate germline workflows under `germline-wgs`; `--library-type` selects the matching configuration, defaulting to `auto-detect` for BAM-based commands.
  - Add `--model-type {large, small}` to `filter_variants germline-wgs` and `train_model germline-wgs`; `large` is the default, and explicit model paths or training iterations take precedence.
- Auto-selected `tumor-normal-wgs` BWA models now use the detected sample type and model-specific score/AF thresholds; explicit CLI values and custom models still take precedence.
- Update the pre-trained cell-line and FFPE BWA models for `tumor-normal-wgs`.
- Turn off base-quality filtering in germline workflows. The `--min-bq` default is now `0` for the `germline`, `germline-multi-sample`, `germline-mitochondria`, and `custom` workflows. The `tumor-normal-wgs`, `tumor-only-wgs`, `germline-tagging`, and `panel-of-normals` workflows keep `--min-bq 23` and now set `--decode-yc consensus` so base types are read from the `YC` tag. The `tumor-only-te` workflow is unchanged.
- Retrain the SBX-Fast germline WGS models with newer sequencing data.
- Base-quality and simplex features remain available but are no longer used by germline models.
- Valid empty BAM/VCF inputs and empty processing results now produce valid header-only outputs with warnings; invalid inputs still fail.
- Change the `--max-indel-size` default for the `tumor-normal-wgs` workflow from `29` to `30`. Indels up to and including 30 bp now pass the indel-size filter.
- BED inputs now accept uncompressed, gzip-compressed, or BGZF-compressed files with a `.gz` extension across Small Variant Caller commands and workflows.
- Reduce peak memory usage during model training.
- Restructured GitBook to multi-page layout.

### Fixed

- Fix multi-allelic germline records with spanning deletions (`*`) so FORMAT fields and `DP` remain valid for downstream tools.
- Record the resolved model and thresholds in the VCF command line for `tumor-normal-wgs`. `filter_variants` now writes the auto-resolved `--model` and per-model `--snv-min-ml-score` and `--indel-min-ml-score` values to the `RocheCommandLine` VCF metadata instead of the pre-resolution default values.
- Preserve BAM features for variants at region boundaries.
- Fix duplicate or stale `POPAF` definitions in output VCF headers.
- Warn when YC decoding encounters an eligible BAM record without a `YC` tag.

### Removed

- Remove the legacy `model-tumor-normal-wgs.txt.gz` and `model-tumor-only-wgs.txt.gz` bundled model aliases; model selection now uses sample-type- and aligner-specific model files.

### Known issues

- Using Giraffe BAMs for germline variant filtering without personalization can lead to reduced performance.

## [1.1.0]

### Added

- Added duplex base-type BAM features computed from the YC tag: `duplex_concordant`, `duplex_simplex`, `duplex_discordant` for ALT alleles and corresponding REF allele features.
- Added VCF FORMAT fields `ADC`, `ADS`, `ADD`, `ADL` (`Number=R`) reporting REF and ALT allelic depths by duplex base type for duplex sequencing workflows.
- Added SHAP value TSV output options in `filter_variants` (`--output-snv-shap-value-tsv`, `--output-indel-shap-value-tsv`, `--output-shap-value-tsv`) for model interpretability. Requires `--target-regions`.
- Added `--max-normal-support`, `--min-tumor-af`, `--max-indel-size`, and `--min-dp-ratio` CLI options to `filter_variants` for the `tumor-normal-wgs` workflow.
- Added `--duplex-lowbq` CLI option in `compute_bam_features` and `filter_variants` to control inclusion of low base-quality reads in duplex depth and allele frequency calculations.
- Added VCF feature `gc_content` for the GC content of the reference sequence in a 200-bp window centered around the variant.
- Added `--aligner` and `--run-type` CLI options to `filter_variants` for `germline` and `germline-multi-sample` workflows to auto-select models based on aligner and run type.
- Added per-model ML score thresholds to the `tumor-normal-wgs` config via a `model_thresholds` map.

### Changed

- Common CLI options are moved from the main command to each workflow-specific subcommand (e.g. `compute_bam_features --warn-as-error germline ...` → `compute_bam_features germline --warn-as-error ...`).
- Standardized TSV and VCF metadata output to use `##RocheCommandLine=<...>` format.
- Improved `filter_variants` performance (~18% wall time reduction) by replacing per-record output assembly with queue-based writer thread and BGZF block-level concatenation.
- Increased default LightGBM training iterations from 3000 to 6000 for the `tumor-normal-wgs` workflow.
- Updated pre-trained cell-line and FFPE models for BWA alignment in the `tumor-normal-wgs` workflow.
- True-positive variants are now excluded from negative training data in `train_model` for the `tumor-normal-wgs` workflow.
- `filter_variants` `tumor-normal-wgs` now rejects `--aligner giraffe` with an error instead of silently falling back to the BWA model. Use `--aligner bwa`, or `--aligner custom` with an explicit `--model`.

### Removed

- Removed `--output-vcf-buffer-size` CLI option from `filter_variants`.
- Removed GATK-prefixed VCF FORMAT fields (`GATK_DP`, `GATK_GT`, `GATK_AD`, `GATK_ALT`) from the germline workflow output.

### Fixed

- Fixed variant density computation producing inconsistent values at region boundaries.
- Fixed incorrect processing of GT fields with three alleles when multi-allelic records are split.
- Fixed `filter_variants` `tumor-normal-wgs` workflow outputting stale FORMAT field values for failed records.
- Fixed integer underflow in `compute_vcf_features` where negative `RPA` INFO values wrapped to large unsigned integers.
- Fixed `compute_vcf_features` POPAF variant mismatch where multi-allelic gnomAD records were not trimmed to their shortest representation.
- Fixed `str` feature always being zero in the `tumor-normal-wgs` workflow.

## [1.0.0]

### Added

- Pre-trained multi-sample models for BWA and Giraffe alignments in the `germline-multi-sample` workflow.
- Support for `tumor-normal-wgs` workflow in `train_model` and `filter_variants`.
- Support for all BAM features in tumor-normal samples.
- Separate `--snv-model` and `--indel-model` options for `filter_variants` in germline workflows.

### Changed

- Replaced `--workflow` CLI option with subcommands for each submodule (`compute_bam_features`, `compute_vcf_features`, `train_model`, `filter_variants`), each with workflow-specific subcommands (e.g. `germline`, `germline-multi-sample`, `tumor-normal-wgs`).
- Renamed numerous CLI options for consistency (see user guide for full mapping).
- Converted boolean CLI flag pairs to alternatives with descriptive values (e.g. `--duplex`/`--no-duplex` → `--sequencing-protocol {duplex, duplex-simplex, umi}`).
- For duplex sequencing protocol, base types are now inferred from the YC tag instead of base quality.
- Reduced peak memory usage for `train_model` by streaming BAM features from disk instead of loading all into memory.
- Deprecated `bam_tumor_af` and `bam_normal_af` features; replaced with `tumor_duplex_af` and `normal_duplex_af`.

### Removed

- Removed `--sample-type` option from `filter_variants`.

### Fixed

- Simplified error messages for invalid input files.
- Fixed splitting and merging of multi-allelic VCF records based on VCF header metadata.
- Fixed parent directory creation for output files.
- Fixed feature normalization to account for sample context in `train_model` and `filter_variants`.

## [0.80.1]

### Added

- `PRED_ML` VCF FORMAT field to report the ML prediction probability score for both `PASS` and `FAIL` records.

### Changed

- `GQ` and `QUAL` fields in output VCF records are now derived from ML prediction probability scores.
- `GT` field is set to `./.` (or `.` for haploid positions) for records with `FAIL` filter instead of `0/1`.
- Updated pre-trained `germline-multi-sample` models for SBX-D and SBX-Fast data.

## [0.80.0]

### Added

- ML-based filtering and re-genotyping of candidate variant calls from GATK HaplotypeCaller for SBX duplex data.
- Workflows: `germline`, `germline-multi-sample`.
- Submodules: `compute_bam_features`, `compute_vcf_features`, `train_model`, `filter_variants`.
- Pre-trained germline models for SBX-D and SBX-Fast chemistries.
- Configurable feature extraction, model training, and variant filtering with JSON config support.

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.0.0...1.1.0
[1.0.0]: https://github.com/Roche-AXELIOS/XOOS/compare/0.80.1...1.0.0
[0.80.1]: https://github.com/Roche-AXELIOS/XOOS/compare/0.80.0...0.80.1
[0.80.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/0.80.0
