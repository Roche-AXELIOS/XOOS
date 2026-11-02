<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Germline Variant Filtering will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- Added fixed two-tier germline classification cutoffs for candidate and population-known variants.
- Added copy-number, loss-of-heterozygosity, and purity-normalized allele-fraction features, together with per-sample segmentation QC metrics.

### Changed

- Renamed training and prediction options for clearer configuration; the release documents the replacement options and compatibility aliases.
- Updated the model feature set and classification defaults. Models trained against an earlier feature layout must be retrained before prediction.
- Restructured the User Guide into a landing README and a multi-page `docs/` set.

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.1.0...1.2.0
