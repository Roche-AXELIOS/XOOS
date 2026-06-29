<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Read Collapser will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0]

### Added

- Duplicate marking via `markdup` subcommand: position-based and UMI-aware clustering of aligned reads, with configurable region padding and BED file input.
- Consensus generation via `consensus` subcommand: error-corrected consensus reads from duplicate clusters, with per-base depth and majority count tags.
- Presets for common workflows: `wgs-duplex`, `wgs-duplex-mrd`, `wgs-duplex-cfdna`, `wgs-simplex`.
- Configurable minimum cluster size, strand-specific cluster size thresholds, and simplex depth filtering.
- Partial read clustering via `--cluster-partials`.
- Merged output mode via `--merge-output` for `markdup`.
- Summary statistics output with detailed read-level metrics (clustering input, unclustered, discarded, and consensus counts).
