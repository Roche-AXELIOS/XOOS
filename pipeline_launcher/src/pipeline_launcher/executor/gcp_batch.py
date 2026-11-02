"""GCP Batch cloud storage and job submission.

Implements the CloudStorage protocol for GCS and GCP Batch. The shared
submission workflow lives in cloud_common.CloudBatchExecutor; this module
provides only the GCP-specific pieces: GCS uploads, job ID sanitization,
and Batch job submission.
"""

from __future__ import annotations

import logging
import mimetypes
import re
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pipeline_launcher.executor.cloud_common import CloudStager

from cloudpathlib import CloudPath
from google.cloud import batch_v1, storage

from pipeline_launcher.config.model import GcpBatchConfig
from pipeline_launcher.executor.base import RunContext, SubmitResult
from pipeline_launcher.executor.formatting import (
    _make_jinja_env,
    _sanitize_label,
)
from pipeline_launcher.executor.path_staging import (
    GcsUri,
    S3Uri,
)

# The gcloud CLI bundled Python — used to override CLOUDSDK_PYTHON when the
# system Python is too old (e.g. Python 3.9 in nextflow/nextflow images).
_CLOUDSDK_BUNDLED_PYTHON = (
    "/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python3"
)

CPU_DRIVER_BOOT_IMAGE = "projects/debian-cloud/global/images/family/debian-12"

# GCP Batch job IDs must match ^[a-z]([a-z0-9-]{0,61}[a-z0-9])?$.
_JOB_ID_INVALID_RE = re.compile(r"[^a-z0-9-]")


def _sanitize_job_id(value: str) -> str:
    """Return a GCP Batch job ID matching ``^[a-z]([a-z0-9-]{0,61}[a-z0-9])?$``."""
    candidate = _JOB_ID_INVALID_RE.sub("-", value.lower())
    candidate = re.sub(r"-+", "-", candidate).strip("-")
    candidate = re.sub(r"^[^a-z]+", "", candidate)
    if not candidate:
        candidate = "job"
    if len(candidate) > 63:
        candidate = candidate[:63].rstrip("-")
    return candidate


# ---------------------------------------------------------------------------
# GcsStorage — CloudStorage implementation for GCP
# ---------------------------------------------------------------------------


class GcsStorage:
    """CloudStorage implementation backed by GCS and GCP Batch."""

    scheme: str = "gs"

    def __init__(self, config: GcpBatchConfig) -> None:
        self._cfg = config.gcp_batch
        self._effective_bucket = self._cfg.gcs_bucket

    @property
    def bucket(self) -> str:
        return self._effective_bucket

    @bucket.setter
    def bucket(self, value: str) -> None:
        self._effective_bucket = value

    def build_uri(self, key: str) -> str:
        return f"gs://{self._effective_bucket}/{key}"

    def upload_file(self, local_path: str, cloud_key: str) -> None:
        client = storage.Client(project=self._cfg.project_id)
        bucket = client.bucket(self._effective_bucket)
        blob = bucket.blob(cloud_key)
        ct, _ = mimetypes.guess_type(local_path)
        blob.upload_from_filename(
            local_path,
            content_type=ct or "application/octet-stream",
        )
        logging.info("Upload successful: gs://%s/%s", self._effective_bucket, cloud_key)

    def upload_directory(self, local_dir: Path, cloud_key_prefix: str) -> int:
        """Upload a directory tree to GCS, skipping hidden files/directories."""
        client = storage.Client(project=self._cfg.project_id)
        bucket = client.bucket(self._effective_bucket)
        count = 0
        for local_file in local_dir.rglob("*"):
            if not local_file.is_file():
                continue
            rel = local_file.relative_to(local_dir)
            if any(part.startswith(".") for part in rel.parts):
                continue
            key = f"{cloud_key_prefix}/{rel.as_posix()}"
            blob = bucket.blob(key)
            ct, _ = mimetypes.guess_type(str(local_file))
            blob.upload_from_filename(
                str(local_file),
                content_type=ct or "application/octet-stream",
            )
            count += 1
        logging.info(
            "Uploaded %d files from %s to gs://%s/%s",
            count,
            local_dir,
            self._effective_bucket,
            cloud_key_prefix,
        )
        return count

    def upload_bytes(self, content: bytes, cloud_key: str, content_type: str) -> None:
        client = storage.Client(project=self._cfg.project_id)
        blob = client.bucket(self._effective_bucket).blob(cloud_key)
        blob.upload_from_string(content, content_type=content_type)

    def render_and_upload_driver(self, template_data: dict, driver_key: str) -> None:
        env = _make_jinja_env()
        template = env.get_template("gcp_batch_driver.jinja")
        bash = template.render(data=template_data)
        self.upload_bytes(bash.encode(), driver_key, "text/x-shellscript")

    def validate_config(self, context: RunContext) -> SubmitResult | None:
        if not self._cfg.gcs_bucket:
            return SubmitResult(
                succeeded=False,
                message="gcp_batch.gcs_bucket must be configured.",
                exit_code=1,
            )
        if not self._cfg.project_id:
            return SubmitResult(
                succeeded=False,
                message="gcp_batch.project_id must be configured.",
                exit_code=1,
            )
        return None

    def make_job_name(self, context: RunContext, unique_id: str) -> str:
        return _sanitize_job_id(
            f"nextflow-driver-{context.analysis_name_unique}-{unique_id}"
        )

    def build_container_command(self, paths: dict[str, str]) -> list[str]:
        # The container downloads the driver from GCS and execs it so
        # Nextflow runs as PID 1 and receives SIGTERM on shutdown.
        return [
            "sh",
            "-c",
            (
                f"export CLOUDSDK_PYTHON={_CLOUDSDK_BUNDLED_PYTHON} && "
                f"DRIVER=$(mktemp) && "
                f"{self._cfg.cli_path} storage cp "
                f"gs://{self._effective_bucket}/{paths['driver_key']} "
                f'"$DRIVER" --quiet && chmod +x "$DRIVER" && '
                f'exec "$DRIVER"'
            ),
        ]

    def get_extra_template_data(self) -> dict[str, str]:
        return {
            "cliPath": self._cfg.cli_path,
            "gcsBucket": self._effective_bucket,
            "cloudsdkPython": _CLOUDSDK_BUNDLED_PYTHON,
        }

    def build_output_url(self, base_key: str) -> str:
        return (
            f"https://console.cloud.google.com/storage/browser/"
            f"{self._effective_bucket}/{base_key}/output"
            f"?project={self._cfg.project_id}"
        )

    def extract_analysis_dir_base(self, analysis_dir: CloudPath) -> str:
        scheme = str(analysis_dir).split("://")[0]
        bucket = analysis_dir.bucket
        return str(analysis_dir)[len(f"{scheme}://{bucket}/") :].rstrip("/")

    def make_stager(self, base_key: str) -> CloudStager:
        from pipeline_launcher.executor.cloud_common import CloudStager

        return CloudStager(
            storage=self,
            native_uri_type=GcsUri,
            foreign_uri_type=S3Uri,
            provider_name="GCP Batch",
            native_scheme_label="gs://",
            foreign_scheme_label="S3",
            base_key=base_key,
        )

    def submit_job(
        self,
        context: RunContext,
        job_name: str,
        container_command: list[str],
    ) -> SubmitResult:
        gcp_cfg = self._cfg
        user = context.username

        runnable = batch_v1.Runnable()
        runnable.container = batch_v1.Runnable.Container()
        runnable.container.image_uri = gcp_cfg.driver_image
        runnable.container.entrypoint = ""
        runnable.container.options = (
            "-v /usr/lib/google-cloud-sdk:/usr/lib/google-cloud-sdk:ro"
        )
        runnable.container.commands = container_command

        task = batch_v1.TaskSpec()
        task.runnables = [runnable]

        group = batch_v1.TaskGroup()
        group.task_count = 1
        group.task_spec = task

        policy = batch_v1.AllocationPolicy.InstancePolicy()
        policy.machine_type = gcp_cfg.driver_machine_type

        boot_disk = batch_v1.AllocationPolicy.Disk()
        boot_disk.image = CPU_DRIVER_BOOT_IMAGE
        boot_disk.size_gb = 10

        policy.boot_disk = boot_disk

        instances = batch_v1.AllocationPolicy.InstancePolicyOrTemplate()
        instances.policy = policy
        instances.install_gpu_drivers = False

        service_account = batch_v1.ServiceAccount(
            email=gcp_cfg.execution_service_account
        )

        allocation_policy = batch_v1.AllocationPolicy()
        allocation_policy.service_account = service_account
        allocation_policy.instances = [instances]

        job = batch_v1.Job()
        job.task_groups = [group]
        job.allocation_policy = allocation_policy
        job.logs_policy.destination = batch_v1.LogsPolicy.Destination.CLOUD_LOGGING
        job.labels = {
            "user": _sanitize_label(user),
            "analysis-name": _sanitize_label(context.analysis_name),
        }

        client = batch_v1.BatchServiceClient()
        create_request = batch_v1.CreateJobRequest(
            parent=f"projects/{gcp_cfg.project_id}/locations/{gcp_cfg.region}",
            job_id=job_name,
            job=job,
        )

        created_job = client.create_job(request=create_request)
        job_name_short = created_job.name.split("/")[-1]
        console_url = (
            f"https://console.cloud.google.com/batch/jobsDetail"
            f"/regions/{gcp_cfg.region}"
            f"/jobs/{job_name_short}"
            f"/logs?project={gcp_cfg.project_id}"
        )

        logging.info("--- Job Submission Successful ---")
        logging.info("Job Name: %s", created_job.name)
        logging.info("Job UID: %s", created_job.uid)
        logging.info("View in GCP Console: %s", console_url)

        return SubmitResult(
            succeeded=True,
            message=f"Submitted GCP Batch job {created_job.name}",
            job_id=created_job.name,  # full resource name; required by poll_job/terminate_job
        )

    def poll_job(self, job_id: str) -> tuple[str, SubmitResult | None]:
        """Check GCP Batch job status; return (status, result).

        Args:
            job_id: Full GCP resource name (projects/{p}/locations/{r}/jobs/{id}).

        Returns:
            Tuple of (status_string, SubmitResult | None). Result is None while
            the job is non-terminal.
        """
        client = batch_v1.BatchServiceClient()
        try:
            job = client.get_job(name=job_id)
        except Exception as exc:
            return "NOT_FOUND", SubmitResult(
                succeeded=False,
                message=f"GCP Batch job {job_id} not found: {exc}",
                job_id=job_id,
                exit_code=1,
            )
        status = job.status.state.name
        if status == "SUCCEEDED":
            return status, SubmitResult(
                succeeded=True,
                message=f"GCP Batch job {job_id} succeeded",
                job_id=job_id,
            )
        if status == "FAILED":
            exit_code = 1
            if job.status.task_groups:
                exec_info = job.status.task_groups[0].task_execution
                if exec_info and exec_info.exit_code:
                    exit_code = exec_info.exit_code
            return status, SubmitResult(
                succeeded=False,
                message=f"GCP Batch job {job_id} failed",
                job_id=job_id,
                exit_code=exit_code,
            )
        return status, None

    def terminate_job(self, job_id: str) -> None:
        """Delete a GCP Batch job (GCP's equivalent of cancel/terminate).

        Args:
            job_id: Full GCP resource name (projects/{p}/locations/{r}/jobs/{id}).
        """
        client = batch_v1.BatchServiceClient()
        client.delete_job(name=job_id)
