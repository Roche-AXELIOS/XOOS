<!-- markdownlint-disable MD024 -->

# Changelog

All notable changes to Pipeline Launcher will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### Added

- Unified CLI tool (`xoos`) for launching Nextflow pipelines across different execution environments.
- Executor backends: `local`, `slurm`, `aws_batch`, and `gcp_batch`.
- `--env-override key=value` CLI option for overriding environment config fields without editing YAML files. Supports dot-notation for nested fields.
- Typer-based CLI with `run` and `rclone` subcommands.
- Parallel rclone copy with file partitioning.
- Cooperative file locking for serializing concurrent pipeline runs.
- Problem report generation (zip of logs and work artifacts) on failure.
- Work directory cleanup with configurable policies.
- Sample status report parsing to determine pipeline success.
- Nextflow config templates and Jinja2 driver templates for authoring environment configs.

<!-- Version comparison links -->
[1.2.0]: https://github.com/Roche-AXELIOS/XOOS/releases/tag/1.2.0
