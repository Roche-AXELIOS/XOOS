"""Slurm executor — renders a driver script and submits via sbatch.

The driver script re-invokes the launcher inside the driver environment
so the actual Nextflow run happens on the allocated compute node with
the correct scratch paths and Singularity cache.
"""

from __future__ import annotations

import json
import logging
import shlex
import subprocess
from pathlib import Path

from pipeline_launcher.config.model import SlurmAttributes, SlurmConfig
from pipeline_launcher.executor.base import RunContext, SubmitResult
from pipeline_launcher.executor.formatting import (
    _make_jinja_env,
    build_driver_command,
    format_shell_command,
)
from pipeline_launcher.paths import stage_dir
from pipeline_launcher.runner.log import rotate


def _build_sbatch_directives(
    job_name: str,
    comment: str,
    attrs: SlurmAttributes,
) -> list[tuple[str, str]]:
    """Return an ordered list of (directive-flag, shell-quoted value) pairs.

    The required ``--job-name`` directive is always first.  Optional
    directives from ``attrs.driver_options`` are appended in iteration
    order; entries whose values are ``None`` or empty are skipped.
    """
    directives: list[tuple[str, str]] = [
        ("--job-name", shlex.quote(job_name)),
    ]
    if comment:
        directives.append(("--comment", shlex.quote(comment)))
    for flag, value in attrs.driver_options.items():
        if value:
            directives.append((f"--{flag}", shlex.quote(value)))
    return directives


def render_slurm_driver_script(
    job_name: str,
    command: list[str],
    comment: str,
    config: SlurmConfig,
    launcher_cwd: str | None = None,
    analysis_dir: str | None = None,
    out_dir: str | None = None,
) -> str:
    """Render a Slurm batch script from the Jinja template."""
    attrs = config.attributes

    data = {
        "sbatchDirectives": _build_sbatch_directives(job_name, comment, attrs),
        "preamble": config.preamble,
        "launcherCwd": shlex.quote(launcher_cwd) if launcher_cwd else None,
        "command": format_shell_command(command),
        "outDir": (
            shlex.quote(out_dir or analysis_dir) if (out_dir or analysis_dir) else None
        ),
    }

    env = _make_jinja_env()
    template = env.get_template("slurm_driver.jinja")
    return template.render(data=data)


def _write_executable_script(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(content)
    path.chmod(0o755)


def _submit_sbatch(script: Path) -> str:
    job_id = subprocess.check_output(
        ["sbatch", "--parsable", str(script)],
        cwd=script.parent,
    )
    return job_id.decode().strip()


class SlurmExecutor:
    """Generate a Slurm driver script and submit it via sbatch."""

    def submit(self, context: RunContext) -> SubmitResult:
        config = context.config
        assert isinstance(config, SlurmConfig)
        assert isinstance(context.out_dir, Path), (
            "Slurm executor requires a local --outdir, got a cloud URI: "
            f"{context.out_dir}"
        )

        if not context.project:
            logging.error("--project must be provided when submitting an HPC job.")
            return SubmitResult(
                succeeded=False,
                message="--project is required for Slurm submission.",
                exit_code=1,
            )

        comment = json.dumps({"project": context.project})

        driver_command = build_driver_command(
            context.analysis_name,
            context.extra_args,
        )

        script_content = render_slurm_driver_script(
            context.analysis_name,
            driver_command,
            comment,
            config,
            launcher_cwd=str(context.launcher_cwd),
            analysis_dir=str(context.analysis_dir),
            out_dir=str(context.out_dir),
        )

        script_path = rotate(stage_dir(context.analysis_dir) / "driver.sh")
        _write_executable_script(script_path, script_content)

        job_id = _submit_sbatch(script_path)
        logging.info(f"Submitted batch job {job_id} with name {context.analysis_name}")

        return SubmitResult(
            succeeded=True,
            message=f"Submitted batch job {job_id}",
            job_id=job_id,
        )
