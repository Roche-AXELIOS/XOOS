<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Pan-Genome Consensus Caller will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0]

### Fixed

- Fixed BED region search returning incorrect results for certain queries.
- Fixed null characters being appended to read names in the clipper output.

## [0.80.1]

### Changed

- Updated documentation.

## [0.80.0]

### Added

- Pangenome-aware consensus correction for duplex reads aligned with Giraffe (`--add-graph-aln`).
- BED file input to restrict processing to non-repetitive regions.
- Multi-threaded processing.
