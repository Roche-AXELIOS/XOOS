"""Environment configuration models.

Each executor type (local, slurm, aws_batch, gcp_batch) has its own model
inheriting from BaseConfig.  The ``executor`` field is used as a
discriminator so ``EnvConfig`` can be parsed directly via
``TypeAdapter[EnvConfig].validate_python(data)``.
"""

from __future__ import annotations

import re
from typing import Annotated, Literal, Union

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    field_validator,
    model_validator,
)

_DRIVER_OPTION_KEY_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9_-]*$")


class _StrictModel(BaseModel):
    """Shared base that rejects unknown fields in all config models."""

    model_config = ConfigDict(extra="forbid")


class SlurmAttributes(_StrictModel):
    """Slurm job scheduler parameters.

    driver_options
        Arbitrary ``#SBATCH`` flags for the driver job that launches
        Nextflow (e.g. ``{"partition": "batch_cpu", "qos": "3d"}``).
        Keys are flag names without the ``--`` prefix; values are the
        flag arguments.  Entries with ``None`` or empty-string values
        are omitted.
    """

    driver_options: dict[str, str | None] = Field(default_factory=dict)

    @field_validator("driver_options")
    @classmethod
    def _validate_driver_options(
        cls, v: dict[str, str | None]
    ) -> dict[str, str | None]:
        for key, value in v.items():
            if not _DRIVER_OPTION_KEY_RE.match(key):
                raise ValueError(
                    f"Invalid driver_options key {key!r}: "
                    f"must match {_DRIVER_OPTION_KEY_RE.pattern}"
                )
            if value is not None and ("\n" in value or "\r" in value):
                raise ValueError(f"driver_options[{key!r}] must not contain newlines")
        return v


class DriverConfig(_StrictModel):
    """HPC driver environment settings for scratch and container caching.

    cache_mode
        Default Singularity cache strategy for this environment:
        ``"shared"`` uses the ``singularity_cache`` path; ``"user"``
        uses a per-user cache under ``scratch_base``.  When unset, the
        CLI default applies.  An explicit ``--singularity-cache`` flag
        always overrides this value.
    """

    scratch_base: str | None = None
    singularity_cache: str | None = None
    cache_mode: Literal["shared", "user"] | None = None

    @model_validator(mode="after")
    def _check_cache_mode_requirements(self) -> "DriverConfig":
        """Fail fast when cache_mode lacks its required path field.

        ``cache_mode: user`` needs ``scratch_base`` (the per-user cache
        lives under it; ``build_nextflow_env`` returns early without it),
        and ``cache_mode: shared`` needs ``singularity_cache``.  Without
        these the strategy would silently not apply at runtime.
        """
        if self.cache_mode == "user" and self.scratch_base is None:
            raise ValueError(
                "driver.cache_mode 'user' requires driver.scratch_base to be set"
            )
        if self.cache_mode == "shared" and self.singularity_cache is None:
            raise ValueError(
                "driver.cache_mode 'shared' requires driver.singularity_cache to be set"
            )
        return self


class AwsBatchAttributes(_StrictModel):
    """AWS Batch job submission parameters."""

    region: str = "us-west-2"
    s3_bucket: str = ""
    driver_queue: str = ""
    cli_path: str = "/mnt/aws-cli/miniconda/bin/aws"
    driver_image: str = "nextflow/nextflow:26.04.0"
    driver_vcpus: int = 2
    driver_memory_mib: int = 3072
    job_role_arn: str = ""
    expected_bucket_owner: str = ""


class GcpBatchAttributes(_StrictModel):
    """GCP Batch job submission parameters."""

    region: str = "us-central1"
    project_id: str = ""
    gcs_bucket: str = ""
    driver_image: str = "nextflow/nextflow:26.04.0"
    execution_service_account: str = ""
    cli_path: str = "/usr/lib/google-cloud-sdk/bin/gcloud"
    driver_machine_type: str = "n2-standard-8"


def _coerce_params(v: dict) -> dict:
    """Coerce YAML booleans/numbers to strings so unquoted values work.

    Accepts a dict whose values may be Python booleans, ints, floats, or
    strings (as produced by a YAML parser) and returns a new dict with all
    values converted to ``str``.  Booleans are lowercased (``True`` →
    ``"true"``) to match Nextflow's convention.  Non-dict input is returned
    unchanged so that Pydantic can surface its own type error.
    """
    if not isinstance(v, dict):
        return v
    return {
        k: str(val).lower() if isinstance(val, bool) else str(val)
        for k, val in v.items()
    }


class BaseConfig(_StrictModel):
    """Fields shared across all executor types."""

    profiles: list[str] = Field(default_factory=list)
    config_files: list[str] = Field(default_factory=list)
    pipeline_params: dict[str, str] = Field(default_factory=dict)

    @field_validator("pipeline_params", mode="before")
    @classmethod
    def _coerce_pipeline_param_values(cls, v: dict) -> dict:
        """Coerce YAML booleans/numbers to strings so unquoted values work."""
        return _coerce_params(v)


class LocalConfig(BaseConfig):
    """Local executor — runs Nextflow directly on the current machine."""

    executor: Literal["local"] = "local"


class SlurmConfig(BaseConfig):
    """Slurm executor — renders an sbatch driver script and submits it."""

    executor: Literal["slurm"] = "slurm"
    attributes: SlurmAttributes = Field(default_factory=SlurmAttributes)
    preamble: str = ""
    driver: DriverConfig = Field(default_factory=DriverConfig)

    @field_validator("preamble", mode="before")
    @classmethod
    def _coerce_preamble(cls, v: str | list[str]) -> str:
        """Accept a list of lines (legacy JSON) or a plain string (YAML)."""
        if isinstance(v, list):
            return "\n".join(v)
        return v


class AwsBatchConfig(BaseConfig):
    """AWS Batch executor — submits a containerized Nextflow driver job."""

    executor: Literal["aws_batch"] = "aws_batch"
    aws_batch: AwsBatchAttributes = Field(default_factory=AwsBatchAttributes)


class GcpBatchConfig(BaseConfig):
    """GCP Batch executor — submits a containerized Nextflow driver job."""

    executor: Literal["gcp_batch"] = "gcp_batch"
    gcp_batch: GcpBatchAttributes = Field(default_factory=GcpBatchAttributes)


class PipelineEnvOverlay(_StrictModel):
    """Pipeline-specific overlay for one env config stem.

    When a pipeline package declares per-environment overrides, each key in its
    ``env`` mapping is validated as a ``PipelineEnvOverlay``.  The overlay is
    merged on top of the base env config at launch time so that pipeline-level
    profiles, config files, and params take precedence over site defaults.
    """

    profiles: list[str] = Field(default_factory=list)
    config_files: list[str] = Field(default_factory=list)
    pipeline_params: dict[str, str] = Field(default_factory=dict)

    @field_validator("pipeline_params", mode="before")
    @classmethod
    def _coerce_pipeline_param_values(cls, v: dict) -> dict:
        """Coerce YAML booleans/numbers to strings so unquoted values work."""
        return _coerce_params(v)


class PipelineConfig(_StrictModel):
    """Top-level pipeline configuration loaded from a pipeline package's config.yaml.

    A pipeline package may ship a *config.yaml* at its root that declares
    default ``pipeline_params`` and optional per-environment overlays under
    the ``env`` key.  The loader validates this file into a ``PipelineConfig``
    instance so that type errors are caught at load time rather than at
    Nextflow invocation.
    """

    pipeline_params: dict[str, str] = Field(default_factory=dict)
    env: dict[str, PipelineEnvOverlay] = Field(default_factory=dict)

    @field_validator("pipeline_params", mode="before")
    @classmethod
    def _coerce_pipeline_param_values(cls, v: dict) -> dict:
        """Coerce YAML booleans/numbers to strings so unquoted values work."""
        return _coerce_params(v)


# Discriminated union — Pydantic selects the right type based on ``executor``.
EnvConfig = Annotated[
    Union[LocalConfig, SlurmConfig, AwsBatchConfig, GcpBatchConfig],
    Field(discriminator="executor"),
]
