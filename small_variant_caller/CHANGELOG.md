<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Small Variant Caller will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
