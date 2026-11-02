"""Build Nextflow CLI commands and environment variables.

Separates command construction from execution so the same logic can be
used by both the local executor (which runs the command directly) and
the slurm executor (which embeds it in a driver script).
"""

from __future__ import annotations

import getpass
import logging
import os
from dataclasses import dataclass
from pathlib import Path

from pipeline_launcher.config.model import DriverConfig
from pipeline_launcher.executor.formatting import build_nf_run_args
from pipeline_launcher.location import LocationPath
from pipeline_launcher.paths import work_dir_target


@dataclass
class PipelineDefinition:
    """Resolved pipeline artifacts."""

    nextflow_script: LocationPath


def create_work_dir(driver: DriverConfig, cwd: Path) -> Path | None:
    """Derive the Nextflow work directory from the driver config.

    Returns None when no scratch is configured, signaling that Nextflow
    should use its default work location.
    """
    if driver.scratch_base is None:
        return None
    scratch = Path(f"{driver.scratch_base}/{getpass.getuser()}")
    scratch_nxf = scratch / ".nextflow"
    scratch_nxf.mkdir(parents=True, exist_ok=True)
    # Mirror the absolute cwd under scratch to keep a stable layout
    # without introducing a leading '/'.
    rel_cwd = Path(*cwd.parts[1:])
    work = scratch_nxf / rel_cwd / "work"
    work.mkdir(parents=True, exist_ok=True)
    return work


def create_work_dir_symlink(out_dir: Path, work: Path) -> None:
    """Create a convenience symlink from out_dir/nextflow/work to the actual work dir.

    Makes the (possibly scratch-hosted) work dir discoverable under the
    output tree for debugging and resuming.

    If the target path already exists we do not overwrite it, but we warn
    when it does not already point at *work*: a real directory (e.g. left by
    a prior local no-scratch run in the same analysis dir) or a symlink to a
    different location would otherwise silently misrepresent where the run's
    work dir actually is, misleading debugging and resume workflows.
    """
    target = work_dir_target(out_dir)
    if not target.exists() and not target.is_symlink():
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(work)
        return

    # Target already exists — only a symlink already pointing at `work` is fine.
    if target.is_symlink() and target.readlink() == work:
        return

    if target.is_symlink():
        logging.warning(
            "Not linking %s to the scratch work dir %s: it is already a symlink "
            "to %s. The displayed work dir may not reflect this run.",
            target,
            work,
            target.readlink(),
        )
    else:
        logging.warning(
            "Not linking %s to the scratch work dir %s: a real path already "
            "exists there (e.g. from a prior local run in this analysis dir). "
            "This run uses %s; the path under the analysis dir is stale.",
            target,
            work,
            work,
        )


def build_nextflow_command(
    definition: PipelineDefinition,
    analysis_name: str,
    config_files: list[Path],
    extra_args: list[str],
    work_dir: Path | None = None,
) -> list[str]:
    """Assemble the full Nextflow CLI argv for local/Slurm execution.

    ``--outdir`` is expected to be already set and canonicalized in
    *extra_args* by the caller (``cli.py`` guarantees this).
    Returns a list starting with ``"nextflow"`` suitable for subprocess use.
    """
    return ["nextflow"] + build_nf_run_args(
        str(definition.nextflow_script),
        analysis_name,
        extra_configs=[str(p) for p in config_files],
        work_dir=str(work_dir) if work_dir is not None else None,
        extra_args=extra_args,
    )


def build_nextflow_env(
    driver: DriverConfig | None,
    singularity_cache: str = "shared",
) -> dict[str, str]:
    """Build environment variables for Nextflow execution.

    Configures JVM sizing, Singularity temp/cache dirs, and the
    NXF_SINGULARITY_CACHEDIR based on the cache strategy.
    """
    cmd_env = dict(os.environ)
    cmd_env["NXF_JVM_ARGS"] = "-Xms1g -Xmx4g"

    if driver is None:
        return cmd_env

    if driver.scratch_base is None:
        return cmd_env

    scratch = Path(f"{driver.scratch_base}/{getpass.getuser()}")

    def _make(suffix: str) -> str:
        p = scratch / suffix
        p.mkdir(parents=True, exist_ok=True)
        return str(p)

    cmd_env["SINGULARITY_TMPDIR"] = _make(".singularity/tmp")
    cmd_env["SINGULARITY_CACHEDIR"] = _make(".singularity/cache")

    if singularity_cache == "user":
        cmd_env["NXF_SINGULARITY_CACHEDIR"] = _make(".nextflow/singularity")
    elif driver.singularity_cache is not None:
        cmd_env["NXF_SINGULARITY_CACHEDIR"] = driver.singularity_cache

    return cmd_env
