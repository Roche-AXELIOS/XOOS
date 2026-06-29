<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Demux will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
