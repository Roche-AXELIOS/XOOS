"""Executor abstraction and factory.

Each executor type implements the same protocol: given a RunContext,
dispatch the pipeline run to the appropriate backend and return a
SubmitResult. The factory reads the config type and returns the
matching executor.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Protocol

from pipeline_launcher.config.model import (
    AwsBatchConfig,
    EnvConfig,
    GcpBatchConfig,
    LocalConfig,
    SlurmConfig,
)
from pipeline_launcher.location import LocationPath
from pipeline_launcher.runner.command import PipelineDefinition


@dataclass
class RunContext:
    """Everything needed to launch a pipeline run."""

    config: EnvConfig
    definition: PipelineDefinition
    analysis_name: str  # Deterministic, hash-free (for all S3/GCS paths)
    out_dir: LocationPath
    launcher_cwd: Path
    analysis_dir: LocationPath | None = None
    analysis_name_unique: str = ""  # Hash-suffixed (for job name/Nextflow only)
    username: str = ""
    env_name: str = ""
    env_raw: str = ""
    extra_args: list[str] = field(default_factory=list)
    lock_path: Path | None = None
    singularity_cache: str = "shared"
    callback: str | None = None
    upload_dst: str | None = None
    rclone_options: str = ""
    work_dir_delete: str = "delete-if-succeeded"
    # Hash-free analysis name base used to build a stable cloud work dir
    # path that survives re-invocations without --analysis-name pinning.
    work_dir_base: str = ""  # Now always set to analysis_name (deterministic)
    project: str | None = None
    pipeline_source_dir: Path | None = None
    # Cloud run sheet: CSV content generated in memory by the CLI.
    # The executor uploads it and adds --input to the Nextflow command.
    run_sheet_content: str | None = None
    # Local samplesheet that needs uploading alongside the run sheet.
    local_samplesheet: Path | None = None
    # When True, cloud executors block until the submitted driver job
    # reaches a terminal state and exit with its status. Off by default so
    # the production launch path stays submit-and-return; only opt-in callers
    # (e.g. the e2e verification driver) set it.
    wait: bool = False
    # Seconds before _poll_until_done gives up waiting for the submitted job.
    # Set from --wait-timeout; default matches the verification driver's 72000s
    # attemptDurationSeconds minus a 600s margin.
    wait_timeout_seconds: int = 71400


@dataclass
class SubmitResult:
    """Outcome of an executor submission."""

    succeeded: bool
    message: str
    job_id: str | None = None
    exit_code: int = 0


class Executor(Protocol):
    """Protocol for pipeline execution backends."""

    def submit(self, context: RunContext) -> SubmitResult: ...


def create_executor(config: EnvConfig) -> Executor:
    """Return the executor matching the config type."""
    # Import here to avoid circular imports and to keep the factory
    # as the single place that knows about all executor types.
    from pipeline_launcher.executor.slurm import SlurmExecutor

    if isinstance(config, LocalConfig):
        from pipeline_launcher.executor.local import LocalDriverExecutor

        return LocalDriverExecutor()
    elif isinstance(config, SlurmConfig):
        return SlurmExecutor()
    elif isinstance(config, AwsBatchConfig):
        from pipeline_launcher.executor.aws_batch import S3Storage
        from pipeline_launcher.executor.cloud_common import CloudBatchExecutor

        return CloudBatchExecutor(S3Storage(config))
    elif isinstance(config, GcpBatchConfig):
        from pipeline_launcher.executor.cloud_common import CloudBatchExecutor
        from pipeline_launcher.executor.gcp_batch import GcsStorage

        return CloudBatchExecutor(GcsStorage(config))
    else:
        raise ValueError(f"No executor for config type {type(config).__name__}")
