<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to MSI Caller will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- Added `tumor-normal` analysis for matched tumor and normal BAMs, which compares repeat-length distributions at microsatellite sites and calculates MSI metrics.
- Added `tumor-only` analysis, which compares tumor repeat-length distributions with a panel-of-normals background model to calculate MSI metrics.
- Added `background_generator` to build reusable JSON background models from panel-of-normals BAMs and threshold-bearing microsatellite sites.

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/compare/1.1.0...1.2.0
