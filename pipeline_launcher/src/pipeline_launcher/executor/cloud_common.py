"""Shared infrastructure for cloud batch executors (AWS Batch, GCP Batch).

Both cloud executors follow the same submission workflow: upload configs
and a driver script to cloud storage, build a Nextflow command, then
submit a batch job that downloads and runs the driver. This module
captures that shared workflow so each provider only implements the
cloud-specific storage and job-submission APIs.

Architecture:
    CloudStorage (Protocol)  — upload/URI operations for a specific provider
    CloudBatchExecutor       — the shared submit() workflow
    CloudStager              — unified passthrough-flag staging (-c, --input, etc.)
    build_cloud_paths()      — cloud-agnostic path layout for a batch run
"""

from __future__ import annotations

import contextlib
import logging
import re
import time
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Protocol

from cloudpathlib import CloudPath

from pipeline_launcher.config.loader import resolve_config_path
from pipeline_launcher.executor.base import RunContext, SubmitResult
from pipeline_launcher.executor.formatting import (
    build_nf_run_args,
    extract_work_dir_from_args,
    format_shell_command,
    render_labels_config,
    resolve_remote_path,
)
from pipeline_launcher.executor.path_staging import (
    DriverFetchFile,
    FlagValue,
    GcsUri,
    LocalFile,
    OpaqueValue,
    S3Uri,
    StagedDataFile,
    stage_allowlisted_args,
)

# Only allow safe characters in directory names used in shell scripts.
_SAFE_DIR_NAME_RE = re.compile(r"^[a-zA-Z0-9._-]+$")

# Fixed path inside the driver container where config files are downloaded.
DRIVER_CONFIG_DIR = "/pipeline/conf"


# ---------------------------------------------------------------------------
# Generic cloud polling loop
# ---------------------------------------------------------------------------


def _poll_until_done(
    *,
    job_id: str,
    poll: Callable[[], tuple[str, SubmitResult | None]],
    terminate: Callable[[], None],
    timeout_seconds: int,
    poll_seconds: int = 30,
) -> SubmitResult:
    """Poll a submitted cloud job until it reaches a terminal state or times out.

    Args:
        job_id: Cloud job identifier, used only in log messages and the timeout result.
        poll: Callable that returns (current_status, result). result is None while
            the job is still running; a SubmitResult when it reaches a terminal state.
        terminate: Best-effort cancel callable, called on KeyboardInterrupt.
            Exceptions raised by terminate are suppressed.
        timeout_seconds: Give up after this many seconds and return a failure result.
        poll_seconds: Seconds to sleep between poll() calls. Use 0 in tests.

    Returns:
        SubmitResult from poll() when terminal, or a failed SubmitResult on timeout.

    Raises:
        KeyboardInterrupt: Re-raised after calling terminate().
    """
    deadline = time.monotonic() + timeout_seconds
    last_status = ""
    logging.info("Waiting for job %s to finish (timeout %ds)", job_id, timeout_seconds)
    try:
        while True:
            status, result = poll()
            elapsed = int(time.monotonic() - (deadline - timeout_seconds))
            if status != last_status:
                logging.info("Job %s: %s (%ds elapsed)", job_id, status, elapsed)
                last_status = status
            if result is not None:
                return result
            if time.monotonic() >= deadline:
                logging.error(
                    "Job %s did not finish within %ds; last status: %s",
                    job_id,
                    timeout_seconds,
                    status,
                )
                return SubmitResult(
                    succeeded=False,
                    message=(
                        f"Job {job_id} did not finish within {timeout_seconds}s "
                        f"(last status: {status})"
                    ),
                    job_id=job_id,
                    exit_code=1,
                )
            time.sleep(poll_seconds)
    except KeyboardInterrupt:
        logging.info("Interrupted — terminating job %s", job_id)
        with contextlib.suppress(Exception):
            terminate()
        raise


# ---------------------------------------------------------------------------
# CloudStorage protocol
# ---------------------------------------------------------------------------


class CloudStorage(Protocol):
    """Provider-specific cloud storage and job submission operations.

    Each cloud provider (S3, GCS) implements this protocol. The
    CloudBatchExecutor delegates all provider-specific work through it.
    """

    scheme: str
    """URI scheme for this provider (``'s3'`` or ``'gs'``)."""

    bucket: str
    """The storage bucket name (mutable for cross-bucket analysis_dir)."""

    def build_uri(self, key: str) -> str:
        """Return the full cloud URI for a storage key."""
        ...

    def upload_file(self, local_path: str, cloud_key: str) -> None:
        """Upload a local file to cloud storage."""
        ...

    def upload_directory(self, local_dir: Path, cloud_key_prefix: str) -> int:
        """Upload a directory tree, skipping hidden files. Returns file count."""
        ...

    def upload_bytes(self, content: bytes, cloud_key: str, content_type: str) -> None:
        """Upload in-memory content to cloud storage."""
        ...

    def render_and_upload_driver(self, template_data: dict, driver_key: str) -> None:
        """Render the provider's Jinja driver template and upload it."""
        ...

    def submit_job(
        self,
        context: RunContext,
        job_name: str,
        container_command: list[str],
    ) -> SubmitResult:
        """Submit the batch job to the cloud provider and return the result."""
        ...

    def poll_job(self, job_id: str) -> tuple[str, SubmitResult | None]:
        """Check the current job status; return (status, result).

        Returns None as the result while the job is still running.
        Returns a SubmitResult when the job has reached a terminal state.
        """
        ...

    def terminate_job(self, job_id: str) -> None:
        """Best-effort cancel. Called on KeyboardInterrupt; exceptions are suppressed by the caller."""
        ...

    def validate_config(self, context: RunContext) -> SubmitResult | None:
        """Validate provider-specific config. Return a failure SubmitResult
        or None if valid."""
        ...

    def make_job_name(self, context: RunContext, unique_id: str) -> str:
        """Generate a provider-appropriate job name."""
        ...

    def build_container_command(self, paths: dict[str, str]) -> list[str]:
        """Build the container entrypoint command that downloads and execs
        the driver script."""
        ...

    def get_extra_template_data(self) -> dict[str, str]:
        """Return provider-specific entries for the driver template data dict."""
        ...

    def build_output_url(self, base_key: str) -> str:
        """Return a provider-specific console URL for the output directory."""
        ...

    def extract_analysis_dir_base(self, analysis_dir: CloudPath) -> str:
        """Extract the key/path portion from a cloud analysis_dir."""
        ...

    def make_stager(self, base_key: str) -> CloudStager:
        """Return a CloudStager configured for this provider."""
        ...


# ---------------------------------------------------------------------------
# Unified path layout
# ---------------------------------------------------------------------------


def build_cloud_paths(
    user: str,
    analysis_name: str,
    bucket: str,
    date: str,
    *,
    scheme: str,
    base: str | None = None,
) -> dict[str, str]:
    """Build the cloud key layout for a batch run.

    *date* should be ``YYYYMMDD``.  When *base* is supplied it overrides
    the auto-computed ``user/date/analysis_name`` layout.
    """
    if base is None:
        base = f"{user}/{date}/{analysis_name}"
    prefix = f"{scheme}://{bucket}/{base}"
    return {
        "base": base,
        "prefix": prefix,
        "work": f"{prefix}/work",
        "output": f"{prefix}/output",
        "report": f"{prefix}/report.html",
        "timeline": f"{prefix}/timeline.html",
        "trace": f"{prefix}/trace.txt",
        "run_log": f"{prefix}/nextflow.log",
        "debug_log": f"{prefix}/.nextflow.log",
        "nextflow_state": f"{prefix}/.nextflow",
        "driver_key": f"{base}/stage/driver.sh",
    }


# ---------------------------------------------------------------------------
# Unified CloudStager
# ---------------------------------------------------------------------------


@dataclass
class CloudStager:
    """Stages passthrough flag values for cloud batch runs.

    Parameterized by the native URI type (the one this cloud accepts)
    and the foreign URI type (the one to reject with an error).
    """

    storage: CloudStorage

    # The URI type native to this cloud (e.g. S3Uri for AWS).
    native_uri_type: type
    # The URI type foreign to this cloud (e.g. GcsUri for AWS).
    foreign_uri_type: type

    # Human-readable labels for error messages.
    provider_name: str  # e.g. "AWS Batch"
    native_scheme_label: str  # e.g. "s3://"
    foreign_scheme_label: str  # e.g. "GCS"

    base_key: str

    def _reject_foreign(self, flag: str, uri: str) -> None:
        raise ValueError(
            f"Flag {flag} value {uri!r} is a {self.foreign_scheme_label} URI "
            f"but this is a {self.provider_name} executor. "
            f"Use a {self.native_scheme_label} URI or a local file path."
        )

    def stage_config(
        self, flag: str, value: FlagValue
    ) -> tuple[str, DriverFetchFile | None]:
        match value:
            case _ if isinstance(value, self.foreign_uri_type):
                uri = value.uri  # type: ignore[union-attr]
                self._reject_foreign(flag, uri)
                # unreachable, but satisfies the type checker
                raise AssertionError  # pragma: no cover
            case S3Uri(uri=uri) | GcsUri(uri=uri):
                filename = uri.rsplit("/", 1)[-1]
                entry = DriverFetchFile(
                    source_uri=uri,
                    container_path=f"{DRIVER_CONFIG_DIR}/{filename}",
                )
                logging.info("Passthrough config %s will be fetched by driver", uri)
                return entry.container_path, entry
            case LocalFile(path=path):
                filename = path.name
                dst_key = f"{self.base_key}/stage/conf/{filename}"
                self.storage.upload_file(str(path), dst_key)
                source_uri = self.storage.build_uri(dst_key)
                entry = DriverFetchFile(
                    source_uri=source_uri,
                    container_path=f"{DRIVER_CONFIG_DIR}/{filename}",
                )
                return entry.container_path, entry
            case OpaqueValue(raw=raw):
                logging.debug(
                    "Config flag %s value %r is not a local file or cloud URI; "
                    "passing through unchanged",
                    flag,
                    raw,
                )
                return raw, None
            case _:
                raise AssertionError(f"Unexpected FlagValue type: {type(value)}")

    def stage_data(
        self, flag: str, value: FlagValue
    ) -> tuple[str, StagedDataFile | None]:
        match value:
            case _ if isinstance(value, self.foreign_uri_type):
                uri = value.uri  # type: ignore[union-attr]
                self._reject_foreign(flag, uri)
                raise AssertionError  # pragma: no cover
            case S3Uri(uri=uri) | GcsUri(uri=uri):
                return uri, StagedDataFile(flag=flag, original_value=uri, cloud_uri=uri)
            case LocalFile(path=path):
                dst_key = f"{self.base_key}/stage/data/{path.name}"
                self.storage.upload_file(str(path), dst_key)
                cloud_uri = self.storage.build_uri(dst_key)
                logging.info("Staged passthrough data %s → %s", str(path), cloud_uri)
                return cloud_uri, StagedDataFile(
                    flag=flag, original_value=str(path), cloud_uri=cloud_uri
                )
            case OpaqueValue(raw=raw):
                logging.debug(
                    "Data flag %s value %r is not a local file or cloud URI; "
                    "passing through unchanged",
                    flag,
                    raw,
                )
                return raw, None
            case _:
                raise AssertionError(f"Unexpected FlagValue type: {type(value)}")


# ---------------------------------------------------------------------------
# CloudRunLayout — groups the resolved paths and bucket for a run
# ---------------------------------------------------------------------------


@dataclass
class CloudRunLayout:
    """Resolved paths, bucket, and job identity for a cloud batch run.

    Groups the values that ``submit()`` computes early and threads
    through every subsequent step, replacing a handful of loose locals.
    """

    effective_bucket: str
    paths: dict[str, str]
    job_name: str
    date: str
    user: str

    @property
    def base_key(self) -> str:
        return self.paths["base"]

    @property
    def driver_key(self) -> str:
        return self.paths["driver_key"]

    @property
    def stage_prefix(self) -> str:
        return f"{self.base_key}/stage"


# ---------------------------------------------------------------------------
# Submit step helpers
# ---------------------------------------------------------------------------


def _resolve_paths_and_bucket(
    context: RunContext,
    storage: CloudStorage,
) -> CloudRunLayout:
    """Compute the cloud path layout and effective bucket for a run.

    When ``analysis_dir`` is already a cloud path, the bucket and base
    key are derived from it so results land alongside the analysis.
    Otherwise the bucket comes from the executor config.
    """
    unique_id = uuid.uuid4().hex[:10]
    job_name = storage.make_job_name(context, unique_id)
    user = context.username
    date = datetime.now(timezone.utc).strftime("%Y%m%d")

    if isinstance(context.analysis_dir, CloudPath):
        effective_bucket = context.analysis_dir.bucket
        dir_base = storage.extract_analysis_dir_base(context.analysis_dir)
        # All subsequent storage.upload_*() calls use storage.bucket
        # internally, so it must reflect the analysis_dir's bucket.
        storage.bucket = effective_bucket
        paths = build_cloud_paths(
            user,
            context.analysis_name,
            effective_bucket,
            date,
            scheme=storage.scheme,
            base=dir_base,
        )
    else:
        effective_bucket = storage.bucket
        paths = build_cloud_paths(
            user,
            context.analysis_name,
            effective_bucket,
            date,
            scheme=storage.scheme,
        )

    return CloudRunLayout(
        effective_bucket=effective_bucket,
        paths=paths,
        job_name=job_name,
        date=date,
        user=user,
    )


def _stage_config_files(
    context: RunContext,
    storage: CloudStorage,
    layout: CloudRunLayout,
) -> list[DriverFetchFile]:
    """Upload resource labels and Nextflow config files to cloud storage.

    Returns the list of DriverFetchFile entries that the driver script
    must download before starting Nextflow.
    """
    # Labels config is generated in-memory and uploaded directly.
    labels_config = render_labels_config(layout.user, context.analysis_name)
    labels_key = f"{layout.base_key}/stage/conf/labels.config"
    storage.upload_bytes(labels_config.encode(), labels_key, "text/plain")
    logging.info("Uploaded labels config to %s", storage.build_uri(labels_key))

    fetch_files: list[DriverFetchFile] = [
        DriverFetchFile(
            source_uri=storage.build_uri(labels_key),
            container_path=f"{DRIVER_CONFIG_DIR}/labels.config",
        ),
    ]

    scheme_prefix = f"{storage.scheme}://"
    foreign_prefix = "gs://" if storage.scheme == "s3" else "s3://"
    for raw_path in context.config.config_files:
        if raw_path.startswith(foreign_prefix):
            raise ValueError(
                f"Config file {raw_path!r} uses the wrong cloud provider. "
                f"Use a {scheme_prefix} URI or a local file path."
            )
        if raw_path.startswith(scheme_prefix):
            # Already in cloud storage — driver fetches directly.
            filename = raw_path.rsplit("/", 1)[-1]
            fetch_files.append(
                DriverFetchFile(
                    source_uri=raw_path,
                    container_path=f"{DRIVER_CONFIG_DIR}/{filename}",
                )
            )
            logging.info("Config %s will be fetched by driver", raw_path)
        else:
            resolved = resolve_config_path(raw_path)
            filename = resolved.name
            dst_key = f"{layout.base_key}/stage/conf/{filename}"
            storage.upload_file(str(resolved), dst_key)
            fetch_files.append(
                DriverFetchFile(
                    source_uri=storage.build_uri(dst_key),
                    container_path=f"{DRIVER_CONFIG_DIR}/{filename}",
                )
            )

    return fetch_files


def _stage_passthrough_args(
    context: RunContext,
    storage: CloudStorage,
    layout: CloudRunLayout,
) -> SubmitResult | list[DriverFetchFile]:
    """Stage the run sheet and allowlisted passthrough args.

    Mutates ``context.extra_args`` in place (adds ``--input`` for the
    run sheet, rewrites staged paths). Returns additional DriverFetchFile
    entries on success, or a failure SubmitResult.
    """
    _stage_run_sheet(context, storage, layout.stage_prefix)

    stager = storage.make_stager(layout.base_key)
    stage_result = stage_allowlisted_args(context.extra_args, stager)
    if isinstance(stage_result, str):
        return SubmitResult(succeeded=False, message=stage_result, exit_code=1)

    context.extra_args = stage_result.rewritten_args
    return list(stage_result.driver_fetch_files)


def _build_nextflow_command(
    context: RunContext,
    storage: CloudStorage,
    layout: CloudRunLayout,
    driver_fetch_files: list[DriverFetchFile],
    pipeline_dir_name: str | None,
) -> str:
    """Build the formatted Nextflow run command for the driver script.

    Resolves the script path, computes the work directory (respecting
    user overrides), and formats the full command string.
    """
    remote_script = _resolve_remote_script(context, pipeline_dir_name)
    driver_config_paths = [f.container_path for f in driver_fetch_files]

    # Respect user-supplied -work-dir; otherwise compute a deterministic
    # path so the same directory is reused across re-invocations.
    passthrough_work_dir = extract_work_dir_from_args(context.extra_args)
    if passthrough_work_dir is not None:
        computed_work_dir = None
    elif isinstance(context.analysis_dir, CloudPath):
        computed_work_dir = layout.paths["work"]
    else:
        scheme_prefix = f"{storage.scheme}://"
        name_base = context.work_dir_base or context.analysis_name
        computed_work_dir = (
            f"{scheme_prefix}{layout.effective_bucket}"
            f"/{layout.user}/{layout.date}/{name_base}/work"
        )

    nf_run = build_nf_run_args(
        script=remote_script,
        analysis_name_unique=context.analysis_name_unique,
        extra_configs=driver_config_paths,
        work_dir=computed_work_dir,
        output_dir=str(context.out_dir),
        extra_args=context.extra_args,
    )
    return format_shell_command(["nextflow"] + nf_run)


def _render_driver(
    storage: CloudStorage,
    layout: CloudRunLayout,
    nf_run_cmd: str,
    driver_fetch_files: list[DriverFetchFile],
    pipeline_source_uri: str | None,
    pipeline_dir_name: str | None,
) -> None:
    """Assemble the template data and render the driver script to cloud storage."""
    template_data: dict[str, Any] = {
        "nextflowCmd": nf_run_cmd,
        "runLogPath": layout.paths["run_log"],
        "debugLogPath": layout.paths["debug_log"],
        "nextflowStatePath": layout.paths["nextflow_state"],
        "driverKey": layout.driver_key,
        "driverFetchFiles": driver_fetch_files,
        "outputUrl": storage.build_output_url(layout.base_key),
    }
    template_data.update(storage.get_extra_template_data())
    if pipeline_source_uri is not None and pipeline_dir_name is not None:
        template_data["pipelineSourceUri"] = pipeline_source_uri
        template_data["pipelineDirName"] = pipeline_dir_name

    storage.render_and_upload_driver(template_data, layout.driver_key)


# ---------------------------------------------------------------------------
# Shared submit workflow
# ---------------------------------------------------------------------------


class CloudBatchExecutor:
    """Executor that delegates cloud-specific operations to a CloudStorage.

    The submit() method orchestrates the submission workflow shared
    between AWS Batch and GCP Batch. Each step is a focused helper;
    provider differences are handled by the injected CloudStorage.
    """

    def __init__(self, storage: CloudStorage) -> None:
        self._storage = storage

    def submit(self, context: RunContext) -> SubmitResult:
        storage = self._storage

        error = storage.validate_config(context)
        if error is not None:
            return error

        layout = _resolve_paths_and_bucket(context, storage)

        driver_fetch_files = _stage_config_files(context, storage, layout)

        pipeline_result = _resolve_pipeline_source(context, storage, layout)
        if isinstance(pipeline_result, SubmitResult):
            return pipeline_result
        pipeline_source_uri, pipeline_dir_name = pipeline_result

        passthrough_result = _stage_passthrough_args(context, storage, layout)
        if isinstance(passthrough_result, SubmitResult):
            return passthrough_result
        driver_fetch_files.extend(passthrough_result)

        nf_run_cmd = _build_nextflow_command(
            context, storage, layout, driver_fetch_files, pipeline_dir_name
        )

        _render_driver(
            storage,
            layout,
            nf_run_cmd,
            driver_fetch_files,
            pipeline_source_uri,
            pipeline_dir_name,
        )

        container_command = storage.build_container_command(layout.paths)
        result = storage.submit_job(
            context,
            layout.job_name,
            container_command,
        )

        # Synchronous mode (--wait): block on the submitted driver job and
        # surface its terminal status. Used by the e2e verification driver,
        # which must fail the verification when the launched pipeline fails.
        if context.wait and result.succeeded and result.job_id is not None:
            return _poll_until_done(
                job_id=result.job_id,
                poll=lambda: storage.poll_job(result.job_id),
                terminate=lambda: storage.terminate_job(result.job_id),
                timeout_seconds=context.wait_timeout_seconds,
            )

        return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_pipeline_source(
    context: RunContext,
    storage: CloudStorage,
    layout: CloudRunLayout,
) -> tuple[str | None, str | None] | SubmitResult:
    """Resolve the pipeline source to a cloud URI or None.

    Three cases:
      1. Local directory  → upload to staging, driver downloads it.
      2. Cloud URI        → already in cloud storage, driver downloads directly.
      3. Neither          → script path baked into the container image.

    Returns (pipeline_source_uri, pipeline_dir_name) on success,
    or a SubmitResult on validation failure.
    """
    script_str = str(context.definition.nextflow_script)
    scheme_prefix = f"{storage.scheme}://"

    if context.pipeline_source_dir is not None:
        pipeline_dir_name = context.pipeline_source_dir.name
        if not _SAFE_DIR_NAME_RE.match(pipeline_dir_name):
            return SubmitResult(
                succeeded=False,
                message=(
                    f"Pipeline directory name {pipeline_dir_name!r} contains "
                    "unsafe characters. Use only alphanumeric, hyphens, "
                    "underscores, and dots."
                ),
                exit_code=1,
            )
        staged_key = f"{layout.base_key}/stage/{pipeline_dir_name}"
        storage.upload_directory(context.pipeline_source_dir, staged_key)
        pipeline_source_uri = storage.build_uri(staged_key)
        return pipeline_source_uri, pipeline_dir_name

    if script_str.startswith(scheme_prefix):
        pipeline_source_uri = script_str.rsplit("/", 1)[0]
        pipeline_dir_name = pipeline_source_uri.rsplit("/", 1)[-1]
        return pipeline_source_uri, pipeline_dir_name

    return None, None


def _resolve_remote_script(
    context: RunContext,
    pipeline_dir_name: str | None,
) -> str:
    """Determine the Nextflow script path for the driver command.

    When a pipeline source directory is involved, paths are rewritten
    to ``${PIPELINE_DIR}/...`` so the driver can resolve them after
    downloading the source to a temp directory.
    """
    script_str = str(context.definition.nextflow_script)
    if pipeline_dir_name is None:
        return script_str

    remote_base = "${PIPELINE_DIR}"
    if context.pipeline_source_dir is not None:
        result = resolve_remote_path(
            context.definition.nextflow_script,
            context.pipeline_source_dir,
            remote_base,
        )
        # nextflow_script is always set, so resolve_remote_path won't
        # return None here.
        assert result is not None
        return result
    return f"{remote_base}/{script_str.rsplit('/', 1)[-1]}"


def _stage_run_sheet(
    context: RunContext,
    storage: CloudStorage,
    stage_prefix: str,
) -> None:
    """Upload the in-memory run sheet CSV (and local samplesheet) to cloud storage.

    The CLI generates a run sheet in memory for cloud runs. This function
    uploads it and appends ``--input <cloud_uri>`` to ``context.extra_args``.
    If a local samplesheet is referenced, it is uploaded first and its
    path is rewritten in the CSV content.
    """
    if context.run_sheet_content is None:
        return

    content = context.run_sheet_content

    if context.local_samplesheet is not None:
        ss = context.local_samplesheet
        ss_key = f"{stage_prefix}/samplesheet/{ss.name}"
        storage.upload_file(str(ss), ss_key)
        ss_uri = storage.build_uri(ss_key)
        content = content.replace(str(ss), ss_uri)

    csv_key = f"{stage_prefix}/run_sheet.csv"
    storage.upload_bytes(content.encode(), csv_key, "text/csv")
    csv_uri = storage.build_uri(csv_key)
    logging.info("Uploaded run sheet to %s", csv_uri)

    context.extra_args.extend(["--input", csv_uri])
