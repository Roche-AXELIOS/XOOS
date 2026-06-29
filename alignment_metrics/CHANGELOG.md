<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Alignment Metrics will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
