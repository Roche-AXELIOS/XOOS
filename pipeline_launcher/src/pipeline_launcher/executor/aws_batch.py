"""AWS Batch cloud storage and job submission.

Implements the CloudStorage protocol for S3 and AWS Batch. The shared
submission workflow lives in cloud_common.CloudBatchExecutor; this module
provides only the AWS-specific pieces: S3 uploads, Batch job definition
management, and job submission.
"""

from __future__ import annotations

import logging
import mimetypes
from pathlib import Path
from typing import Any

import boto3
import boto3.s3.transfer
from cloudpathlib import CloudPath

from pipeline_launcher.config.model import AwsBatchConfig
from pipeline_launcher.executor.base import RunContext, SubmitResult
from pipeline_launcher.executor.cloud_common import CloudStager
from pipeline_launcher.executor.formatting import (
    _make_jinja_env,
    _sanitize_label,
)
from pipeline_launcher.executor.path_staging import (
    GcsUri,
    S3Uri,
)

# Job definition name prefix used when creating definitions dynamically.
_JOB_DEF_PREFIX = "xoos-nf-driver"


def _guess_content_type(path: str) -> str:
    """Guess MIME type from a file path, defaulting to binary."""
    suffix = Path(path).suffix.lower()
    if suffix in (".config", ".log"):
        return "text/plain"
    ct, _ = mimetypes.guess_type(path)
    return ct or "application/octet-stream"


# ---------------------------------------------------------------------------
# AWS Batch job definition management
# ---------------------------------------------------------------------------


def _find_matching_job_definition(
    batch_client: Any,
    jd_name: str,
    desired_image: str,
    cli_mount: str,
    job_role_arn: str,
) -> str | None:
    """Return the ARN of an existing active job definition whose immutable
    properties match.

    Only ``image``, the AWS-CLI volume ``sourcePath``, and ``jobRoleArn``
    are checked — these require a new revision when changed. Properties
    overridable at submit time (e.g. resourceRequirements) are excluded.
    """
    response = batch_client.describe_job_definitions(
        jobDefinitionName=jd_name, status="ACTIVE"
    )
    definitions = sorted(
        response.get("jobDefinitions", []),
        key=lambda d: d.get("revision", 0),
        reverse=True,
    )
    for jd in definitions:
        props = jd.get("containerProperties", {})
        if props.get("image") != desired_image:
            continue
        volumes = props.get("volumes", [])
        cli_volume = next((v for v in volumes if v.get("name") == "aws-cli"), None)
        if cli_volume is None:
            continue
        if cli_volume.get("host", {}).get("sourcePath") != cli_mount:
            continue
        existing_role = props.get("jobRoleArn", "")
        if existing_role != job_role_arn:
            continue
        return jd["jobDefinitionArn"]
    return None


def _ensure_job_definition(
    batch_client: Any,
    batch_cfg: Any,
) -> str:
    """Return the ARN of a Batch job definition matching the current config.

    Reuses an existing active revision when immutable properties match,
    avoiding runaway revision growth on repeated launches.
    """
    jd_name = _JOB_DEF_PREFIX

    # The AWS CLI is bind-mounted from the host. The mount covers the
    # install root (two levels above the binary) so shared libraries
    # are available inside the container.
    cli_mount = str(Path(batch_cfg.cli_path).parent.parent)

    existing_arn = _find_matching_job_definition(
        batch_client, jd_name, batch_cfg.driver_image, cli_mount, batch_cfg.job_role_arn
    )
    if existing_arn is not None:
        logging.info("Reusing job definition %s", existing_arn)
        return existing_arn

    container_props: dict = {
        "image": batch_cfg.driver_image,
        "resourceRequirements": [
            {"type": "VCPU", "value": str(batch_cfg.driver_vcpus)},
            {"type": "MEMORY", "value": str(batch_cfg.driver_memory_mib)},
        ],
        "volumes": [
            {"host": {"sourcePath": cli_mount}, "name": "aws-cli"},
        ],
        "mountPoints": [
            {
                "containerPath": cli_mount,
                "readOnly": True,
                "sourceVolume": "aws-cli",
            },
        ],
    }
    if batch_cfg.job_role_arn:
        container_props["jobRoleArn"] = batch_cfg.job_role_arn

    updated = batch_client.register_job_definition(
        jobDefinitionName=jd_name,
        type="container",
        containerProperties=container_props,
    )
    jd_arn = updated["jobDefinitionArn"]
    logging.info("Registered job definition %s", jd_arn)
    return jd_arn


def _exit_code_from_job(job: dict) -> int:
    """Extract the container exit code from a FAILED Batch job description.

    Falls back to 1 when no attempts or exit code are recorded.

    Args:
        job: A job description dict from ``describe_jobs``.

    Returns:
        The exit code from the last attempt's container, or 1 if unavailable.
    """
    attempts = job.get("attempts") or []
    if attempts:
        return attempts[-1].get("container", {}).get("exitCode", 1) or 1
    return 1


# ---------------------------------------------------------------------------
# S3Storage — CloudStorage implementation for AWS
# ---------------------------------------------------------------------------


class S3Storage:
    """CloudStorage implementation backed by S3 and AWS Batch."""

    scheme: str = "s3"

    def __init__(self, config: AwsBatchConfig) -> None:
        self._cfg = config.aws_batch
        self._effective_bucket = self._cfg.s3_bucket
        self._output_prefix = ""

    @property
    def bucket(self) -> str:
        return self._effective_bucket

    @bucket.setter
    def bucket(self, value: str) -> None:
        self._effective_bucket = value

    def build_uri(self, key: str) -> str:
        return f"s3://{self._effective_bucket}/{key}"

    def upload_file(self, local_path: str, cloud_key: str) -> None:
        s3 = boto3.client("s3")
        ct = _guess_content_type(local_path)
        extra_args: dict = {"ContentType": ct}
        owner = self._cfg.expected_bucket_owner
        if owner:
            extra_args["ExpectedBucketOwner"] = owner
        s3.upload_file(
            local_path, self._effective_bucket, cloud_key, ExtraArgs=extra_args
        )
        logging.info("Upload successful: s3://%s/%s", self._effective_bucket, cloud_key)

    def upload_directory(self, local_dir: Path, cloud_key_prefix: str) -> int:
        """Upload a directory tree to S3, skipping hidden files/directories."""
        s3 = boto3.client("s3")
        transfer = boto3.s3.transfer.S3Transfer(s3)
        owner = self._cfg.expected_bucket_owner
        count = 0
        for local_file in local_dir.rglob("*"):
            if not local_file.is_file():
                continue
            relative = local_file.relative_to(local_dir)
            if any(part.startswith(".") for part in relative.parts):
                continue
            s3_key = f"{cloud_key_prefix}/{relative.as_posix()}"
            ct = _guess_content_type(str(local_file))
            extra_args: dict = {"ContentType": ct}
            if owner:
                extra_args["ExpectedBucketOwner"] = owner
            transfer.upload_file(
                str(local_file),
                self._effective_bucket,
                s3_key,
                extra_args=extra_args,
            )
            count += 1
        logging.info(
            "Uploaded %d files from %s to s3://%s/%s",
            count,
            local_dir,
            self._effective_bucket,
            cloud_key_prefix,
        )
        return count

    def upload_bytes(self, content: bytes, cloud_key: str, content_type: str) -> None:
        s3 = boto3.client("s3")
        kwargs: dict = {
            "Bucket": self._effective_bucket,
            "Key": cloud_key,
            "Body": content,
            "ContentType": content_type,
        }
        owner = self._cfg.expected_bucket_owner
        if owner:
            kwargs["ExpectedBucketOwner"] = owner
        s3.put_object(**kwargs)

    def render_and_upload_driver(self, template_data: dict, driver_key: str) -> None:
        env = _make_jinja_env()
        template = env.get_template("aws_batch_driver.jinja")
        bash = template.render(data=template_data)
        self.upload_bytes(bash.encode(), driver_key, "text/x-shellscript")

    def validate_config(self, context: RunContext) -> SubmitResult | None:
        if not self._cfg.s3_bucket:
            return SubmitResult(
                succeeded=False,
                message="aws_batch.s3_bucket must be configured.",
                exit_code=1,
            )
        return None

    def make_job_name(self, context: RunContext, unique_id: str) -> str:
        return f"driver-analysis-{context.analysis_name_unique}-{unique_id}"

    def build_container_command(self, paths: dict[str, str]) -> list[str]:
        # Stash the base prefix for the S3 console URL in submit_job.
        self._output_prefix = paths["base"]

        # The container downloads the driver from S3 and execs it so
        # Nextflow runs as PID 1 and receives SIGTERM on shutdown.
        return [
            "sh",
            "-c",
            (
                f"DRIVER=$(mktemp) && "
                f"{self._cfg.cli_path} s3 cp "
                f"s3://{self._effective_bucket}/{paths['driver_key']} "
                f'"$DRIVER" --quiet && chmod +x "$DRIVER" && '
                f'exec "$DRIVER"'
            ),
        ]

    def get_extra_template_data(self) -> dict[str, str]:
        return {
            "cliDir": str(Path(self._cfg.cli_path).parent),
            "s3Bucket": self._effective_bucket,
        }

    def build_output_url(self, base_key: str) -> str:
        region = self._cfg.region
        return (
            f"https://{region}.console.aws.amazon.com/s3/buckets/"
            f"{self._effective_bucket}?region={region}"
            f"&prefix={base_key}/output/&showversions=false"
        )

    def extract_analysis_dir_base(self, analysis_dir: CloudPath) -> str:
        return analysis_dir.key.rstrip("/")

    def make_stager(self, base_key: str) -> CloudStager:
        from pipeline_launcher.executor.cloud_common import CloudStager

        return CloudStager(
            storage=self,
            native_uri_type=S3Uri,
            foreign_uri_type=GcsUri,
            provider_name="AWS Batch",
            native_scheme_label="s3://",
            foreign_scheme_label="GCS",
            base_key=base_key,
        )

    def submit_job(
        self,
        context: RunContext,
        job_name: str,
        container_command: list[str],
    ) -> SubmitResult:
        batch_cfg = self._cfg
        batch_client = boto3.client("batch", region_name=batch_cfg.region)

        job_def = _ensure_job_definition(batch_client, batch_cfg)

        container_overrides = {
            "command": container_command,
            "resourceRequirements": [
                {"type": "VCPU", "value": str(batch_cfg.driver_vcpus)},
                {"type": "MEMORY", "value": str(batch_cfg.driver_memory_mib)},
            ],
        }

        user = context.username
        response = batch_client.submit_job(
            jobName=job_name,
            jobQueue=batch_cfg.driver_queue,
            jobDefinition=job_def,
            containerOverrides=container_overrides,
            tags={
                "user": _sanitize_label(user),
                "analysis-name": _sanitize_label(context.analysis_name),
            },
        )

        job_id = response["jobId"]
        job_url = (
            f"https://console.aws.amazon.com/batch/home"
            f"?region={batch_cfg.region}#/jobs/ec2/detail/{job_id}"
        )
        output_url = (
            f"https://{batch_cfg.region}.console.aws.amazon.com/s3/buckets/"
            f"{self._effective_bucket}?region={batch_cfg.region}"
            f"&prefix={self._output_prefix}/&showversions=false"
        )

        logging.info("--- Job Submission Successful ---")
        logging.info("Job Name: %s", response["jobName"])
        logging.info("Job ID: %s", job_id)
        logging.info("View in AWS Batch Console: %s", job_url)
        logging.info("View results: %s", output_url)

        return SubmitResult(
            succeeded=True,
            message=f"Submitted AWS Batch job {response['jobName']}",
            job_id=job_id,
        )

    def poll_job(self, job_id: str) -> tuple[str, SubmitResult | None]:
        """Check AWS Batch job status; return (status, result).

        Returns None as the result while the job is non-terminal.
        Maps SUCCEEDED and FAILED to terminal SubmitResults.

        Args:
            job_id: The AWS Batch job ID.

        Returns:
            Tuple of (status_string, SubmitResult | None).
        """
        batch_client = boto3.client("batch", region_name=self._cfg.region)
        jobs = batch_client.describe_jobs(jobs=[job_id]).get("jobs", [])
        if not jobs:
            return "NOT_FOUND", SubmitResult(
                succeeded=False,
                message=f"AWS Batch job {job_id} not found",
                job_id=job_id,
                exit_code=1,
            )
        job = jobs[0]
        status = job.get("status", "UNKNOWN")
        if status == "SUCCEEDED":
            return status, SubmitResult(
                succeeded=True,
                message=f"AWS Batch job {job_id} succeeded",
                job_id=job_id,
            )
        if status == "FAILED":
            reason = job.get("statusReason", "unknown")
            return status, SubmitResult(
                succeeded=False,
                message=f"AWS Batch job {job_id} failed: {reason}",
                job_id=job_id,
                exit_code=_exit_code_from_job(job),
            )
        return status, None

    def terminate_job(self, job_id: str) -> None:
        """Terminate an AWS Batch job. Works for all states including RUNNING.

        Uses terminate_job (not cancel_job) because cancel_job silently no-ops for
        RUNNING jobs.

        Args:
            job_id: The AWS Batch job ID.
        """
        batch_client = boto3.client("batch", region_name=self._cfg.region)
        batch_client.terminate_job(jobId=job_id, reason="Parent process terminated")
