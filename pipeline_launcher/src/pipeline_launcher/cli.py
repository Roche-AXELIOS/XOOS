"""CLI entry point for the pipeline launcher.

Provides two subcommands:
  - run:   Launch a Nextflow pipeline (dispatches to local, slurm,
           aws_batch, or gcp_batch)
  - rclone: Parallel rclone copy using file partitioning

The executor type is determined by the environment config, not by the
subcommand. The user always runs `xoos run --env <name>` regardless
of whether the backend is local, Slurm, AWS Batch, or GCP Batch.
"""

from __future__ import annotations

import logging
import re
import signal
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple, Optional

import click
import typer
from cloudpathlib import CloudPath

from pipeline_launcher.location import LocationPath, to_location

if TYPE_CHECKING:
    from pipeline_launcher.config.model import BaseConfig

_KEY_SEGMENT_RE = re.compile(r"^[a-zA-Z_]\w*$")


class KeyValueParam(click.ParamType):
    """Click parameter type that parses ``key=value`` strings.

    The key must be one or more valid Python identifier segments separated
    by dots (e.g. ``driver.singularity_cache``).  Everything after the
    first ``=`` is the value (may be empty or contain additional ``=``).
    """

    name = "KEY=VALUE"

    def convert(
        self,
        value: str,
        param: click.Parameter | None,
        ctx: click.Context | None,
    ) -> tuple[str, str]:
        key, sep, val = value.partition("=")
        if not sep:
            self.fail(
                f"Expected KEY=VALUE format, got {value!r}. "
                f"Example: driver.singularity_cache=/my/path",
                param,
                ctx,
            )
        for segment in key.split("."):
            if not _KEY_SEGMENT_RE.match(segment):
                self.fail(
                    f"Invalid key {key!r}: each segment must be a valid "
                    f"Python identifier (got {segment!r}).",
                    param,
                    ctx,
                )
        return (key, val)


KEY_VALUE = KeyValueParam()

app = typer.Typer(
    name="xoos",
    help="Launch Nextflow pipelines across local, Slurm, AWS Batch, and GCP Batch environments.",
    add_completion=False,
)


def _setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)-8s - %(message)s",
    )


def _setup_signal_handler() -> None:
    """Convert SIGINT and SIGTERM into KeyboardInterrupt for a single cancellation path.

    Handling SIGTERM is required for driver mode: when a Slurm job is cancelled via
    ``scancel``, Slurm sends SIGTERM (not SIGINT) to the driver process.  Without this
    handler Python would terminate immediately and the graceful Nextflow shutdown in
    ``run_with_log`` (``process.terminate()`` / ``process.kill()``) would never execute.
    """

    def handler(signum: int, frame: object) -> None:
        raise KeyboardInterrupt

    signal.signal(signal.SIGINT, handler)
    signal.signal(signal.SIGTERM, handler)


def _passthrough_param_keys(passthrough_args: list[str]) -> set[str]:  # noqa: E302
    """Extract parameter names from passthrough args (e.g. --outdir -> outdir)."""
    keys: set[str] = set()
    for arg in passthrough_args:
        if arg.startswith("--"):
            # Handle both --key value and --key=value forms.
            key = arg.lstrip("-").split("=", 1)[0]
            keys.add(key)
    return keys


def _extract_passthrough_profile(passthrough_args: list[str]) -> str | None:
    """Extract and consume ``-profile`` from passthrough args.

    If the user passed ``-profile foo,bar`` in the passthrough args,
    return ``"foo,bar"`` and remove both the flag and its value from
    the list so they are not forwarded twice.  Returns ``None`` when
    no ``-profile`` is present.
    """
    for i, arg in enumerate(passthrough_args):
        if arg == "-profile" and i + 1 < len(passthrough_args):
            value = passthrough_args[i + 1]
            del passthrough_args[i : i + 2]
            return value
    return None


def _build_profile_args(
    env_config: "BaseConfig",
    passthrough: list[str],
) -> tuple[list[str], list[str]]:
    """Merge env and user-supplied profiles into a single ``-profile`` arg.

    Consumes ``-profile`` from *passthrough* (mutating the list in-place)
    and deduplicates the combined list while preserving order.

    Returns ``(profile_args, passthrough)`` where ``passthrough`` has the
    ``-profile`` flag and its value removed.
    """
    user_profile = _extract_passthrough_profile(passthrough)
    all_profiles = list(env_config.profiles) + (
        user_profile.split(",") if user_profile else []
    )
    seen: set[str] = set()
    profiles: list[str] = []
    for p in all_profiles:
        if p not in seen:
            seen.add(p)
            profiles.append(p)
    if profiles:
        return ["-profile", ",".join(profiles)], passthrough
    return [], passthrough


def _build_pipeline_param_args(
    env_config: "BaseConfig",
    override_keys: set[str],
) -> list[str]:
    """Build ``--key value`` args for pipeline params not overridden by passthrough.

    Parameters from ``passthrough_args`` take precedence: any key present
    in ``override_keys`` is omitted so Nextflow sees only the user-supplied
    value.
    """
    args: list[str] = []
    for k, v in env_config.pipeline_params.items():
        if k not in override_keys:
            args.extend([f"--{k}", str(v)])
    return args


def _build_resources_base_args(
    resources_base: LocationPath | str | None,
    config_files: list[str] | None,
) -> list[str]:
    """Inject ``--resources_base`` and route ``resources.config``.

    For cloud executors (*config_files* is not ``None``) the
    ``resources.config`` path is appended to *config_files* so the
    executor can stage it into the driver container.  For local/Slurm
    runs it is passed as a ``-c`` flag and its existence is validated.

    Raises ``SystemExit(1)`` when the local ``resources.config`` is missing.
    """
    if resources_base is None:
        return []

    resources_base_loc = (
        to_location(resources_base)
        if isinstance(resources_base, str)
        else resources_base
    )
    resources_base_str = str(resources_base_loc)
    resources_config = f"{resources_base_str}/resources.config"

    args: list[str] = ["--resources_base", resources_base_str]

    if config_files is not None:
        config_files.append(resources_config)
    else:
        if (
            not isinstance(resources_base_loc, CloudPath)
            and not Path(resources_config).is_file()
        ):
            logging.error(
                "resources.config not found at %s. "
                "Ensure --resources-base points to a directory containing resources.config.",
                resources_config,
            )
            raise SystemExit(1)
        args.extend(["-c", resources_config])

    return args


def _build_bundled_config_args(
    config_files: list[str] | None,
    *,
    disable_provenance_report: bool = False,
    disable_sample_report: bool = False,
    enable_pb_dynamic_scaling: bool = False,
    disable_per_sample_error_ignore: bool = False,
) -> list[str]:
    """Inject bundled Nextflow config files and load required plugins.

    ``nf_prov.config`` and ``nf_report.config`` are included unless
    explicitly disabled.  ``pb_dynamic_scaling.config`` is only included
    when explicitly enabled via ``--enable pb-dynamic-scaling``.

    ``fail_on_ignore.config`` is included by default (sets
    ``workflow.failOnIgnore = true``) and suppressed when
    ``--disable per-sample-error-ignore`` switches the error strategy to ``finish``.

    Nextflow ignores ``plugins { }`` blocks in configs supplied via ``-c``,
    so plugins are loaded explicitly via the ``-plugins`` flag instead.

    For cloud executors (*config_files* is not ``None``) the bundled
    config names are appended to *config_files* so the executor uploads
    them.  For local/Slurm runs they are passed as ``-c`` flags.
    In both cases the ``-plugins`` flag is returned so it flows through
    to the Nextflow command.
    """
    from pipeline_launcher.config.loader import bundled_nfconfig

    bundled: list[str] = []
    # Always load nf-schema for config validation in the pipeline
    plugins: list[str] = ["nf-schema@2.7.2"]
    if not disable_provenance_report:
        bundled.append("@nextflow_config/nf_prov.config")
        plugins.append("nf-prov@1.7.0")
    if not disable_sample_report:
        bundled.append("@nextflow_config/nf_report.config")
        plugins.append("nf-report@1.1.1")
    if enable_pb_dynamic_scaling:
        bundled.append("@nextflow_config/pb_dynamic_scaling.config")
    # fail_on_ignore.config is appended last so it takes precedence
    # over any pipeline-level failOnIgnore setting (launcher-wins).
    if not disable_per_sample_error_ignore:
        bundled.append("@nextflow_config/fail_on_ignore.config")

    plugin_args = ["-plugins", ",".join(plugins)] if plugins else []

    if config_files is not None:
        config_files.extend(bundled)
        return plugin_args

    args: list[str] = []
    for b in bundled:
        args.extend(["-c", str(bundled_nfconfig(b.removeprefix("@nextflow_config/")))])
    return args + plugin_args


def _build_error_strategy_args(
    disable_per_sample_error_ignore: bool,
    override_keys: set[str],
) -> list[str]:
    """Build the ``--error_strategy_fallback`` Nextflow param.

    Only injected when ``disable_per_sample_error_ignore`` is ``True``,
    switching the fallback from the env config default (``ignore``) to
    ``finish``.  In the default case no param is injected — the env
    configs use ``params.getOrDefault('error_strategy_fallback', 'ignore')``.

    Skipped if the user already supplied ``--error_strategy_fallback``
    in passthrough args or pipeline_params.
    """
    if not disable_per_sample_error_ignore:
        return []
    if "error_strategy_fallback" in override_keys:
        return []
    return ["--error_strategy_fallback", "finish"]


def _build_extra_args(
    env_config: "BaseConfig",
    passthrough_args: list[str] | None = None,
    *,
    resources_base: LocationPath | str | None = None,
    disable_provenance_report: bool = False,
    disable_sample_report: bool = False,
    enable_pb_dynamic_scaling: bool = False,
    disable_per_sample_error_ignore: bool = False,
    config_files: list[str] | None = None,
) -> tuple[list[str], list[str]]:
    """Assemble Nextflow extra arguments from config and CLI flags.

    Returns a tuple of ``(extra_args, filtered_passthrough)`` where
    ``filtered_passthrough`` is the passthrough list with consumed
    arguments (like ``-profile``) removed.

    Parameters from ``passthrough_args`` take precedence: any key that
    appears in passthrough args is omitted from pipeline_params so
    Nextflow sees only the user-supplied value.

    If the user passed ``-profile`` in passthrough args, those profiles
    are appended to the environment profiles and emitted as a single
    ``-profile`` argument.

    Scheduler attributes for Slurm (partition, qos, account) live under
    ``attributes.driver_options`` and apply to the driver sbatch job
    only.  To control the account charged to individual Nextflow task
    jobs, set ``slurm_account`` in ``pipeline_params``; it is forwarded
    as a ``--slurm_account`` Nextflow param and consumed by the env
    config's ``process.clusterOptions``.
    """
    passthrough = list(passthrough_args or [])
    passthrough_keys = _passthrough_param_keys(passthrough)
    # Merge passthrough and pipeline_params keys so that launcher-injected
    # params (e.g. error_strategy_fallback) are suppressed when either
    # source already provides the same key.
    override_keys = passthrough_keys | set(env_config.pipeline_params)

    profile_args, passthrough = _build_profile_args(env_config, passthrough)
    param_args = _build_pipeline_param_args(env_config, passthrough_keys)
    resources_args = _build_resources_base_args(resources_base, config_files)
    error_strategy_args = _build_error_strategy_args(
        disable_per_sample_error_ignore, override_keys
    )
    bundled_args = _build_bundled_config_args(
        config_files,
        disable_provenance_report=disable_provenance_report,
        disable_sample_report=disable_sample_report,
        enable_pb_dynamic_scaling=enable_pb_dynamic_scaling,
        disable_per_sample_error_ignore=disable_per_sample_error_ignore,
    )

    extra_args = (
        profile_args + param_args + resources_args + error_strategy_args + bundled_args
    )
    return extra_args, passthrough


@dataclass
class RunParams:
    """Grouped CLI parameters for a pipeline run."""

    env: str
    pipeline_script: LocationPath
    output_dir: LocationPath
    launcher_cwd: Path
    username: str
    analysis_dir: LocationPath | None = None
    resources_base: LocationPath | None = None
    upload_dst: str | None = None
    file_lock: Path | None = None
    rclone_options: str = ""
    name: str | None = None
    singularity_cache: str | None = None
    callback: str | None = None
    work_dir_delete: str = "delete-if-succeeded"
    project: str | None = None
    analysis_name: str | None = None
    driver_mode: bool = False
    disable: list[str] = field(default_factory=list)
    enable: list[str] = field(default_factory=list)
    passthrough_args: list[str] = field(default_factory=list)
    run_sheet_content: str | None = None
    local_samplesheet: Path | None = None
    env_overrides: list[tuple[str, str]] = field(default_factory=list)
    pipeline: str | None = None
    pipeline_explicit: bool = False
    wait: bool = False
    wait_timeout_seconds: int = 71400


class _RunSheetResult(NamedTuple):
    """Result of run sheet generation."""

    passthrough: list[str]
    # CSV content for cloud upload (None for local runs).
    run_sheet_content: str | None
    # Local samplesheet path that the executor must upload (None when
    # the samplesheet is already a cloud URI or for local runs).
    local_samplesheet: Path | None


@dataclass
class RunSheetParams:
    """Optional run-sheet generation parameters for ``xoos run``.

    When any field is provided, all five must be supplied together.
    Mutually exclusive with ``--input`` in the passthrough args.
    """

    run_name: str | None = None
    run_type: str | None = None
    file_type: str | None = None
    run_dir: str | None = None
    samplesheet: str | None = None
    workflow_config: str | None = None
    resources_base: str | None = None


def _generate_run_sheet_if_needed(
    passthrough: list[str],
    run_sheet: RunSheetParams,
    *,
    is_cloud: bool = False,
) -> _RunSheetResult:
    """Generate a run sheet if run-sheet params are provided.

    For local/Slurm runs the CSV is written to the output directory and
    ``--input <path>`` is prepended to the passthrough args.

    For cloud runs the CSV content is returned in-memory via
    ``run_sheet_content`` so the executor can upload it directly.  If
    the samplesheet is a local file, ``local_samplesheet`` is set so
    the executor can upload it alongside the run sheet.
    """
    run_sheet_args = [
        run_sheet.run_type,
        run_sheet.file_type,
        run_sheet.run_dir,
        run_sheet.samplesheet,
    ]
    provided = [a for a in run_sheet_args if a is not None]

    if not provided:
        # Check if --input is in passthrough args; if not, throw an error since we need either --input or run sheet args to proceed.
        if "--input" not in passthrough and not any(
            arg.startswith("--input=") for arg in passthrough
        ):
            logging.error(
                "Either --input must be provided in passthrough args or all of "
                "--run-name, --run-type, --file-type, --run-dir, and "
                "--samplesheet must be provided to generate a run sheet."
            )
            raise SystemExit(1)
        return _RunSheetResult(passthrough, None, None)

    # All five must be provided together.
    if len(provided) != len(run_sheet_args):
        missing = []
        for flag, val in [
            ("--run-name", run_sheet.run_name),
            ("--run-type", run_sheet.run_type),
            ("--file-type", run_sheet.file_type),
            ("--run-dir", run_sheet.run_dir),
            ("--samplesheet", run_sheet.samplesheet),
        ]:
            if val is None:
                missing.append(flag)
        logging.error(
            "When generating a run sheet, all five flags are required. Missing: %s",
            ", ".join(missing),
        )
        raise SystemExit(1)

    # --input must not also appear in passthrough args.
    if "--input" in passthrough:
        logging.error(
            "--input in passthrough args conflicts with run sheet "
            "generation flags (--run-name, etc.). Use one or the other."
        )
        raise SystemExit(1)

    from pipeline_launcher.samplesheet import (
        RunSheetRow,
        generate_run_sheet,
        generate_run_sheet_content,
    )

    row = RunSheetRow(
        run_name=run_sheet.run_name,  # type: ignore[arg-type]
        run_type=run_sheet.run_type,  # type: ignore[arg-type]
        file_type=run_sheet.file_type,  # type: ignore[arg-type]
        run_dir=to_location(run_sheet.run_dir),  # type: ignore[arg-type]
        samplesheet=to_location(run_sheet.samplesheet),  # type: ignore[arg-type]
    )

    if is_cloud:
        content = generate_run_sheet_content(row)
        # If the samplesheet is a local file, the executor needs to
        # upload it.  The run sheet content already references the
        # resolved local path; the executor will rewrite it after
        # uploading.
        local_ss: Path | None = None
        if not isinstance(to_location(str(run_sheet.samplesheet)), CloudPath):
            local_ss = Path(str(run_sheet.samplesheet)).resolve()
        return _RunSheetResult(passthrough, content, local_ss)

    out_dir = _extract_outdir(passthrough)
    assert isinstance(out_dir, Path), (
        "local run sheet generation requires a local --outdir"
    )
    sheet_path = generate_run_sheet(row, out_dir)
    return _RunSheetResult(["--input", str(sheet_path)] + passthrough, None, None)


def _extract_outdir(args: list[str]) -> LocationPath | None:
    """Extract --outdir from passthrough args, or return ``None`` if absent.

    Returns a ``CloudPath`` for cloud URIs (preserving ``s3://`` scheme) or
    a resolved absolute ``Path`` for local directories.

    Raises ``SystemExit(1)`` when ``--outdir`` is present but has no value
    (e.g. it is the last argument).
    """
    for i, arg in enumerate(args):
        if arg == "--outdir":
            if i + 1 >= len(args):
                logging.error("--outdir requires a value")
                raise SystemExit(1)
            return to_location(args[i + 1])
        if arg.startswith("--outdir="):
            return to_location(arg.split("=", 1)[1])
    return None


def _canonicalize_outdir(args: list[str]) -> None:
    """Convert --outdir path in passthrough args to an absolute path.

    Cloud URIs (s3://, gs://) are left unchanged.
    """
    for i, arg in enumerate(args):
        if arg == "--outdir" and i + 1 < len(args):
            loc = to_location(args[i + 1])
            if not isinstance(loc, CloudPath):
                args[i + 1] = str(loc)
        elif arg.startswith("--outdir="):
            val = arg.split("=", 1)[1]
            loc = to_location(val)
            if not isinstance(loc, CloudPath):
                args[i] = f"--outdir={loc}"


def _resolve_username(explicit: str | None) -> str:
    """Return *explicit* if set, otherwise detect from the OS."""
    if explicit:
        return explicit
    import getpass
    import os

    try:
        return getpass.getuser()
    except (KeyError, OSError):
        return os.environ.get("USER", f"uid-{os.getuid()}")


def _validate_outdir_locality(out_dir: LocationPath, *, is_cloud: bool) -> None:
    """Ensure --outdir locality matches the executor type."""
    outdir_is_cloud = isinstance(out_dir, CloudPath)
    if is_cloud and not outdir_is_cloud:
        logging.error(
            "--outdir must be a cloud URI (s3://, gs://, az://) "
            "when using a cloud executor, got: %s",
            out_dir,
        )
        raise SystemExit(1)
    if not is_cloud and outdir_is_cloud:
        logging.error(
            "--outdir must be a local path when using a local/Slurm executor, got: %s",
            out_dir,
        )
        raise SystemExit(1)


_DISABLE_VALID = frozenset(
    {
        "provenance-report",
        "sample-report",
        "per-sample-error-ignore",
    }
)


def _parse_disable_flags(disable: list[str]) -> tuple[bool, bool, bool]:
    """Parse ``--disable`` values into individual boolean flags.

    Accepted values: ``provenance-report``, ``sample-report``,
    ``per-sample-error-ignore``.  Unknown values cause a logged error
    and ``SystemExit(1)``.

    Returns ``(disable_provenance_report, disable_sample_report,
    disable_per_sample_error_ignore)``.
    """
    unknown = [v for v in disable if v not in _DISABLE_VALID]
    if unknown:
        logging.error(
            "Unknown --disable value(s): %s. Valid values: %s",
            ", ".join(unknown),
            ", ".join(sorted(_DISABLE_VALID)),
        )
        raise SystemExit(1)
    return (
        "provenance-report" in disable,
        "sample-report" in disable,
        "per-sample-error-ignore" in disable,
    )


_ENABLE_VALID = frozenset(
    {
        "pb-dynamic-scaling",
    }
)


def _parse_enable_flags(enable: list[str]) -> tuple[bool]:
    """Parse ``--enable`` values into individual boolean flags.

    Accepted values: ``pb-dynamic-scaling``.  Unknown values cause a
    logged error and ``SystemExit(1)``.

    Returns ``(enable_pb_dynamic_scaling,)``.
    """
    unknown = [v for v in enable if v not in _ENABLE_VALID]
    if unknown:
        logging.error(
            "Unknown --enable value(s): %s. Valid values: %s",
            ", ".join(unknown),
            ", ".join(sorted(_ENABLE_VALID)),
        )
        raise SystemExit(1)
    return ("pb-dynamic-scaling" in enable,)


def _load_and_override_config(
    env: str, overrides: list[tuple[str, str]]
) -> tuple["BaseConfig", Path]:
    """Load env config, apply CLI overrides, and return (config, path)."""
    from pipeline_launcher.config.loader import (
        apply_env_overrides,
        load_env_config,
        resolve_env_config_path,
    )

    env_config_path = resolve_env_config_path(env)
    env_config = load_env_config(env_config_path)
    if overrides:
        try:
            apply_env_overrides(env_config, overrides)
        except ValueError as exc:
            logging.error("--env-override error: %s", exc)
            raise SystemExit(1) from exc
    return env_config, env_config_path


def _resolve_singularity_cache(cli_value: str | None, env_config: object) -> str:
    """Resolve the effective Singularity cache strategy.

    Precedence: explicit ``--singularity-cache`` flag, then the env
    config's ``driver.cache_mode``, then the ``"shared"`` default.
    """
    if cli_value is not None:
        if cli_value not in ("shared", "user"):
            logging.error(
                "Invalid --singularity-cache %r: must be 'shared' or 'user'.",
                cli_value,
            )
            raise SystemExit(1)
        return cli_value

    driver = getattr(env_config, "driver", None)
    cache_mode = getattr(driver, "cache_mode", None) if driver else None
    return cache_mode or "shared"


def _resolve_and_apply_pipeline_config(
    env_config: "BaseConfig",
    pipeline: str | None,
    pipeline_source_dir: "Path | None",
    env_name: str,
    *,
    is_explicit: bool = True,
) -> None:
    """Resolve the pipeline config and apply it to env_config in-place.

    No-op when no pipeline config can be resolved (generic / undetectable).

    Raises:
        ValueError: When the pipeline config YAML is malformed or invalid.

    Args:
        env_config: Env config to mutate (profiles, config_files, pipeline_params).
        pipeline: Effective pipeline name or path (already resolved by caller).
        pipeline_source_dir: Directory containing the pipeline's nextflow.config.
        env_name: Stem of the env config path, used for overlay lookup.
        is_explicit: Forwarded to ``resolve_pipeline_config_path``; controls
            whether an unresolvable name aborts or falls back to generic.
    """
    from pipeline_launcher.config.loader import (
        apply_pipeline_config,
        load_pipeline_config,
        resolve_pipeline_config_path,
    )

    pipeline_config_path = resolve_pipeline_config_path(
        pipeline, pipeline_source_dir, is_explicit=is_explicit
    )
    if pipeline_config_path is None:
        return
    pipeline_cfg = load_pipeline_config(pipeline_config_path)
    apply_pipeline_config(env_config, pipeline_cfg, env_name)


def _resolve_effective_pipeline(
    pipeline: str | None,
    pipeline_script: str,
) -> str | None:
    """Return the effective pipeline name, auto-detecting from manifest.name if needed.

    When ``pipeline`` is ``None``, reads ``nextflow.config`` from the pipeline
    script's parent directory (local or cloud) and extracts ``manifest.name``.
    Returns the resolved name, or ``None`` when detection fails.

    Args:
        pipeline: The explicit ``--pipeline`` value, or ``None`` for auto-detect.
        pipeline_script: Path or cloud URI to the Nextflow pipeline script.

    Returns:
        The resolved pipeline name string, or ``None``.
    """
    if pipeline is not None:
        return pipeline

    from pipeline_launcher.config.loader import detect_pipeline_name

    loc = to_location(pipeline_script)
    if isinstance(loc, CloudPath):
        detected = detect_pipeline_name(loc.parent)
    else:
        script_path = Path(str(loc))
        parent = script_path.parent
        if not script_path.is_file() or not (parent / "nextflow.config").exists():
            return None
        detected = detect_pipeline_name(parent)

    if detected:
        logging.debug("Auto-detected pipeline: %s", detected)
    return detected


def _resolve_pipeline_script(
    raw_script: LocationPath,
) -> tuple[LocationPath, "Path | None"]:
    """Resolve a pipeline script to a usable reference and detect its source directory.

    Cloud URIs are returned as-is with no source directory.  Local paths are
    resolved to absolute and validated; the parent directory is returned as the
    source directory when it contains a sibling ``nextflow.config``.

    Args:
        raw_script: The pipeline script path or cloud URI from ``RunParams``.

    Returns:
        ``(script_ref, pipeline_source_dir)`` — ``pipeline_source_dir`` is
        the parent directory of local scripts with a sibling ``nextflow.config``,
        or ``None`` for cloud scripts and scripts without one.

    Raises:
        SystemExit: When a local script path does not exist.
    """
    if isinstance(raw_script, CloudPath):
        logging.debug("Using cloud pipeline script: %s", raw_script)
        return raw_script, None
    script_path = Path(str(raw_script)).resolve()
    if not script_path.is_file():
        logging.error("Pipeline script not found: %s", script_path)
        raise SystemExit(1)
    parent = script_path.parent
    if (parent / "nextflow.config").exists():
        logging.debug("Detected local pipeline directory at %s", parent)
        return script_path, parent
    return script_path, None


def _is_bundled_pipeline_name(name: str | None) -> bool:
    """Return True when *name* identifies a bundled pipeline (not generic, not a file path).

    Args:
        name: Candidate pipeline name or ``None``.

    Returns:
        ``True`` when the name is a non-generic string that does not resolve to
        an existing file path on disk.
    """
    return bool(name and name != "generic" and not Path(name).exists())


def _execute_run(params: RunParams) -> None:
    """Resolve config, build context, and dispatch to the executor."""
    from pipeline_launcher.executor.base import RunContext, create_executor
    from pipeline_launcher.naming import (
        generate_analysis_name,
        generate_analysis_name_unique,
    )
    from pipeline_launcher.runner.command import PipelineDefinition

    env_config, env_config_path = _load_and_override_config(
        params.env, params.env_overrides
    )

    singularity_cache = _resolve_singularity_cache(params.singularity_cache, env_config)

    # For cloud executors, bundled configs are appended to config_files
    # so the executor uploads them alongside env configs. For
    # local/Slurm they are passed as -c flags directly.
    from pipeline_launcher.config.model import AwsBatchConfig, GcpBatchConfig

    is_cloud = isinstance(env_config, (AwsBatchConfig, GcpBatchConfig))

    (
        disable_provenance_report,
        disable_sample_report,
        disable_per_sample_error_ignore,
    ) = _parse_disable_flags(params.disable)

    (enable_pb_dynamic_scaling,) = _parse_enable_flags(params.enable)

    # Resolve the pipeline script.  Cloud URIs are passed through as-is;
    # local paths are resolved to absolute and checked for a parent
    # directory containing nextflow.config (which cloud executors upload).
    script_ref, pipeline_source_dir = _resolve_pipeline_script(params.pipeline_script)

    # Apply pipeline-specific config overlays before building extra args so
    # that pipeline_params, profiles, and config_files are fully composed
    # when _build_extra_args reads env_config.
    try:
        _resolve_and_apply_pipeline_config(
            env_config,
            params.pipeline,
            pipeline_source_dir,
            env_config_path.stem,
            is_explicit=params.pipeline_explicit,
        )
    except ValueError as exc:
        logging.error("Pipeline config error: %s", exc)
        raise SystemExit(1) from exc

    if params.driver_mode:
        # The passthrough args were fully assembled by the first invocation
        # (profiles merged, bundled configs added, -plugins injected, etc.).
        # Re-running _build_extra_args here would duplicate every
        # launcher-injected arg (-plugins, -c config files, --resources_base),
        # causing Nextflow to receive conflicting flags and fail to apply the
        # Slurm task-executor profile correctly.
        extra_args = list(params.passthrough_args)
    else:
        extra_args, remaining_passthrough = _build_extra_args(
            env_config,
            params.passthrough_args,
            resources_base=params.resources_base,
            disable_provenance_report=disable_provenance_report,
            disable_sample_report=disable_sample_report,
            enable_pb_dynamic_scaling=enable_pb_dynamic_scaling,
            disable_per_sample_error_ignore=disable_per_sample_error_ignore,
            config_files=env_config.config_files if is_cloud else None,
        )
        extra_args.extend(remaining_passthrough)

    analysis_name = generate_analysis_name(
        name_override=params.name,
        analysis_name_override=params.analysis_name,
    )
    analysis_name_unique = generate_analysis_name_unique(
        name_override=params.name,
        analysis_name_override=params.analysis_name,
    )

    definition = PipelineDefinition(
        nextflow_script=script_ref,
    )

    context = RunContext(
        config=env_config,
        definition=definition,
        analysis_name=analysis_name,
        analysis_name_unique=analysis_name_unique,
        out_dir=params.output_dir,
        launcher_cwd=params.launcher_cwd,
        username=params.username,
        env_name=env_config_path.stem if params.env else "",
        env_raw=params.env or "",
        extra_args=extra_args,
        lock_path=params.file_lock,
        singularity_cache=singularity_cache,
        callback=params.callback,
        upload_dst=params.upload_dst,
        rclone_options=params.rclone_options,
        work_dir_delete=params.work_dir_delete,
        work_dir_base=analysis_name,
        project=params.project,
        pipeline_source_dir=pipeline_source_dir,
        run_sheet_content=params.run_sheet_content,
        local_samplesheet=params.local_samplesheet,
        analysis_dir=params.analysis_dir,
        wait=params.wait,
        wait_timeout_seconds=params.wait_timeout_seconds,
    )

    # In driver mode (re-invoked inside a Slurm job), always use the
    # local executor so Nextflow runs directly instead of submitting
    # another Slurm job.
    if params.driver_mode:
        from pipeline_launcher.executor.local import LocalExecutor

        executor = LocalExecutor()
    else:
        executor = create_executor(env_config)
    result = executor.submit(context)

    if not result.succeeded:
        raise SystemExit(1)


class _WorkflowConfigResult(NamedTuple):
    """Values derived from a workflow configuration YAML."""

    run_name: str | None
    run_type: str | None
    file_type: str | None
    run_dir: str | None
    samplesheet: str | None


def _apply_workflow_config(  # noqa: PLR0913
    workflow_config: str,
    *,
    run_name: str | None,
    run_type: str | None,
    file_type: str | None,
    run_dir: str | None,
    samplesheet: str | None,
    resources_base: str | None,
    launcher_cwd: Path,
    passthrough: list[str],
) -> _WorkflowConfigResult:
    """Parse a workflow config YAML and derive run parameters.

    Explicit CLI values (non-None) are preserved; only None fields are
    filled from the YAML.  *passthrough* is mutated in place to inject
    derived analysis params and the ``--workflow_config`` provenance arg.

    Raises ``ValueError`` on invalid YAML or derivation errors.
    """
    from pipeline_launcher.workflow_config import (
        derive_file_type,
        derive_pipeline_params,
        derive_run_dir,
        derive_run_type,
        generate_samplesheet_csv,
        load_workflow_config,
        warn_reference_integrity,
    )

    wf = load_workflow_config(to_location(workflow_config))

    # Derive defaults — explicit CLI flags win.
    if run_name is None:
        run_name = wf.workflow_order_name
    if run_type is None:
        run_type = derive_run_type(wf)
    if file_type is None:
        file_type = derive_file_type(wf)
    if run_dir is None:
        run_dir = str(derive_run_dir(wf))

    # Generate samplesheet from Samples[] if not provided.
    # Always write to a local staging directory — for cloud executors
    # the existing _RunSheetResult.local_samplesheet mechanism handles
    # uploading.
    if samplesheet is None and wf.samples:
        stage_dir = launcher_cwd / ".xoos_stage"
        samplesheet = str(generate_samplesheet_csv(wf, stage_dir))

    # Inject analysis params as low-priority passthrough args.
    derived_params = derive_pipeline_params(wf)
    existing_keys = _passthrough_param_keys(passthrough)
    for k, v in derived_params.items():
        if k not in existing_keys:
            passthrough.extend([f"--{k}", v])

    # Reference genome existence check.
    warn_reference_integrity(
        wf,
        to_location(resources_base) if resources_base is not None else None,
    )

    # Forward the YAML path to Nextflow for provenance.
    existing_keys = _passthrough_param_keys(passthrough)
    if "workflow_config" not in existing_keys:
        passthrough.extend(["--workflow_config", workflow_config])

    return _WorkflowConfigResult(
        run_name=run_name,
        run_type=run_type,
        file_type=file_type,
        run_dir=run_dir,
        samplesheet=samplesheet,
    )


@app.command(
    context_settings={"allow_extra_args": True, "allow_interspersed_args": False},
)
def run(  # noqa: PLR0913
    ctx: typer.Context,  # NOSONAR — Typer requires one parameter per CLI option
    env: str = typer.Option(
        ...,
        help=(
            "Environment name or path to an env config JSON. "
            "Looked up as env/{name}.json in the package; "
            "falls back to a direct file path."
        ),
    ),
    pipeline_script: str = typer.Option(
        ...,
        help=(
            "Path to the Nextflow pipeline script (.nf), or a cloud "
            "URI (s3://, gs://). Local paths inside a directory with "
            "nextflow.config are uploaded by cloud executors."
        ),
    ),
    resources_base: str = typer.Option(
        ...,
        help=(
            "Base path for pipeline resources (local directory or "
            "cloud URI like s3://bucket/resources). Must contain a "
            "resources.config file."
        ),
    ),
    # Run-sheet generation — optional alternative to --input in passthrough.
    run_name: Optional[str] = typer.Option(
        None,
        help=(
            "Run name for auto-generated run sheet. "
            "When set, --run-type, --file-type, --run-dir, and "
            "--samplesheet are also required and --input must not "
            "appear in passthrough args."
        ),
    ),
    run_type: Optional[str] = typer.Option(
        None,
        help="Run type (e.g. SBX-D, SBX-FAST, YSU).",
    ),
    file_type: Optional[str] = typer.Option(
        None,
        help="Input file type (e.g. basecall_rdb, basecall_fastq, demultiplexed_fastq, bam, cram).",
    ),
    run_dir: Optional[str] = typer.Option(
        None,
        help="Path or cloud URI (s3://, gs://) to the directory containing run data.",
    ),
    samplesheet: Optional[str] = typer.Option(
        None,
        help="Path or cloud URI to the per-run samplesheet CSV (sample_name, sample_sid, ...).",
    ),
    upload_dst: Optional[str] = typer.Option(
        None,
        help="Destination for uploading results via rclone.",
    ),
    file_lock: Optional[Path] = typer.Option(
        None,
        help="Path to a lock file for preventing concurrent executions.",
    ),
    rclone_options: str = typer.Option(
        "",
        help=(
            "Options to pass to rclone (quoted string), "
            "e.g. '--transfers=20 --copy-links'."
        ),
    ),
    name: Optional[str] = typer.Option(
        None,
        help="Base name for analysis (used for job name and Nextflow -name).",
    ),
    analysis_dir: Optional[str] = typer.Option(
        None,
        help=(
            "Required. Root directory (or cloud URI) for all pipeline artifacts "
            "(stage, nextflow dir, work dir, logs). "
            "--outdir defaults to {analysis-dir}/output when not provided."
        ),
    ),
    singularity_cache: Optional[str] = typer.Option(
        None,
        help=(
            "Singularity cache strategy: 'shared' or 'user'. Overrides the "
            "env config's driver.cache_mode. Defaults to driver.cache_mode "
            "when set, otherwise 'shared'."
        ),
    ),
    callback: Optional[str] = typer.Option(
        None,
        help="Command to invoke with START/COMPLETE/FAIL around the Nextflow run.",
    ),
    work_dir_delete: str = typer.Option(
        "delete-if-succeeded",
        help="Work dir cleanup: delete-if-succeeded, delete-always, delete-never.",
    ),
    project: Optional[str] = typer.Option(
        None,
        help="Project name for Slurm job comment. Required for HPC submission.",
    ),
    disable: list[str] = typer.Option(
        [],
        help=(
            "Features to disable. Accepted values (repeatable): "
            "provenance-report, sample-report, "
            "per-sample-error-ignore (switches default errorStrategy from 'ignore' "
            "to 'finish' and suppresses workflow.failOnIgnore)."
        ),
    ),
    enable: list[str] = typer.Option(
        [],
        help=(
            "Features to enable. Accepted values (repeatable): "
            "pb-dynamic-scaling (inject input-size-aware memory scaling "
            "for Parabricks processes)."
        ),
    ),
    username: Optional[str] = typer.Option(
        None,
        help=(
            "Override the username used for cloud work directory paths "
            "and resource labels. Defaults to the current OS user."
        ),
    ),
    # Hidden args used by the Slurm driver to preserve submitter context.
    launcher_cwd: Path = typer.Option(
        Path.cwd().absolute(),
        hidden=True,
    ),
    analysis_name: Optional[str] = typer.Option(
        None,
        hidden=True,
    ),
    driver_mode: bool = typer.Option(
        False,
        hidden=True,
        help="Internal flag set by the Slurm driver to force local execution.",
    ),
    # Annotated as list[str] because Typer cannot introspect
    # list[tuple[str, str]].  At runtime, click_type=KEY_VALUE makes each
    # element a (key, value) tuple — see KeyValueParam.convert().
    env_override: Optional[list[str]] = typer.Option(
        None,
        click_type=KEY_VALUE,
        help=(
            "Override env config fields. Repeatable. Format: key=value "
            "with dot-notation for nested fields "
            "(e.g. --env-override driver.singularity_cache=/my/path)."
        ),
    ),
    workflow_config: Optional[str] = typer.Option(
        None,
        help=(
            "Path to a workflow configuration YAML. "
            "Derives --run-name, --run-type, --file-type, --run-dir, "
            "and --samplesheet from the YAML. Explicit CLI flags "
            "override derived values."
        ),
    ),
    pipeline: Optional[str] = typer.Option(
        None,
        help=(
            "Pipeline name or path to a pipeline config YAML. "
            "Bundled pipelines: 'xoosnf'. "
            "Use 'generic' to skip pipeline-specific config. "
            "Defaults to auto-detect from the pipeline script's manifest.name."
        ),
    ),
    wait: bool = typer.Option(
        False,
        "--wait",
        help=(
            "Block until the submitted cloud driver job finishes and exit "
            "with its status, instead of returning immediately after "
            "submission. Cloud executors only (AWS Batch). Intended for "
            "synchronous callers such as the e2e verification driver."
        ),
    ),
    wait_timeout: int = typer.Option(
        71400,
        "--wait-timeout",
        help=(
            "Seconds before --wait gives up polling the submitted cloud job. "
            "Defaults to 71400 (20 h minus a 10-minute margin). "
            "Only meaningful when --wait is also set."
        ),
    ),
) -> None:
    """Launch a Nextflow pipeline on the configured environment.

    Extra arguments after ``--`` are forwarded directly to Nextflow.

    --analysis-dir sets the root directory (or cloud URI) for all pipeline artifacts (stage, .nextflow, work, logs, output). If not provided, it is auto-derived per environment.

    --outdir is a Nextflow-only parameter. If specified, it is passed through as-is. If not, it is set to {analysis_dir}/output.

    -profile can be included as a passthrough argument to specify additional Nextflow profiles defined by the pipeline or additional config files supplied via -c.

    Example (passthrough --input)::

        xoos run --env /path/to/env.yaml -- \
            --analysis-dir /analysis/root --input sheet.csv -profile germline_wgs

    Example (generated run sheet)::

        xoos run --env /path/to/env.yaml \
            --run-name my_run --run-type SBX-D \
            --file-type basecall_fastq --run-dir /data/run_001 \
            --samplesheet /data/samples.csv -- \
            --analysis-dir /analysis/root -profile germline_wgs

    Example (workflow config — derives run sheet params from YAML)::

        xoos run --env /path/to/env.yaml \
            --workflow-config /data/run_001/workflow-configuration.yaml -- \
            --analysis-dir /analysis/root -profile germline_wgs
    """
    _setup_logging()
    _setup_signal_handler()

    try:
        from pipeline_launcher.config.loader import load_env_config as _load
        from pipeline_launcher.config.loader import resolve_env_config_path
        from pipeline_launcher.config.model import AwsBatchConfig, GcpBatchConfig

        _env_cfg = _load(resolve_env_config_path(env))
        is_cloud = isinstance(_env_cfg, (AwsBatchConfig, GcpBatchConfig))

        passthrough = list(ctx.args)

        # --analysis-dir is required for all executors (local, Slurm, cloud).
        # Every executor uses it to locate the stage dir, nextflow dir, and
        # work dir — there is no sensible default to guess.
        if analysis_dir is None:
            logging.error(
                "--analysis-dir is required. "
                "Provide a local path (local/Slurm) or cloud URI (cloud executors), "
                "e.g. --analysis-dir /runs/my-run  or  --analysis-dir s3://bucket/user/run"
            )
            raise SystemExit(1)
        analysis_dir_val: LocationPath = to_location(analysis_dir)

        # --outdir defaults to {analysis_dir}/output; the user can override
        # it by passing --outdir explicitly after --.
        # This must happen before run sheet generation because local run
        # sheet generation writes the CSV into the output directory.
        if _extract_outdir(passthrough) is None:
            passthrough.extend(
                ["--outdir", f"{str(analysis_dir_val).rstrip('/')}/output"]
            )

        # Resolve the effective pipeline name before run-sheet generation so
        # both plugin loading and _execute_run's config application use the
        # same name — including when auto-detection is triggered.
        _effective_pipeline = _resolve_effective_pipeline(pipeline, pipeline_script)

        # --- Plugin discovery: use the resolved effective pipeline name ---
        _plugin = None
        if _is_bundled_pipeline_name(_effective_pipeline):
            from pipeline_launcher.pipeline import load_run_sheet_plugin

            _plugin = load_run_sheet_plugin(_effective_pipeline)

        rs_params = RunSheetParams(
            run_name=run_name,
            run_type=run_type,
            file_type=file_type,
            run_dir=run_dir,
            samplesheet=samplesheet,
            workflow_config=workflow_config,
            resources_base=resources_base,
        )

        if _plugin is not None:
            # The plugin handles both workflow config derivation and run sheet
            # generation internally.
            try:
                rs_result = _plugin.generate(
                    params=rs_params,
                    launcher_cwd=Path(launcher_cwd),
                    is_cloud=is_cloud,
                    passthrough=passthrough,
                )
            except ValueError as exc:
                logging.error("Plugin run sheet error: %s", exc)
                raise SystemExit(1) from exc
        else:
            # --- Workflow config: derive defaults from workflow YAML ---
            if workflow_config is not None:
                try:
                    wf_result = _apply_workflow_config(
                        workflow_config,
                        run_name=run_name,
                        run_type=run_type,
                        file_type=file_type,
                        run_dir=run_dir,
                        samplesheet=samplesheet,
                        resources_base=resources_base,
                        launcher_cwd=Path(launcher_cwd),
                        passthrough=passthrough,
                    )
                    run_name = wf_result.run_name
                    run_type = wf_result.run_type
                    file_type = wf_result.file_type
                    run_dir = wf_result.run_dir
                    samplesheet = wf_result.samplesheet
                except ValueError as exc:
                    logging.error("--workflow-config error: %s", exc)
                    raise SystemExit(1) from exc

            rs_result = _generate_run_sheet_if_needed(
                passthrough,
                RunSheetParams(
                    run_name=run_name,
                    run_type=run_type,
                    file_type=file_type,
                    run_dir=run_dir,
                    samplesheet=samplesheet,
                ),
                is_cloud=is_cloud,
            )

        passthrough = rs_result.passthrough

        out_dir = _extract_outdir(passthrough)
        assert out_dir is not None  # guaranteed: injected above if absent
        _validate_outdir_locality(out_dir, is_cloud=is_cloud)
        _canonicalize_outdir(passthrough)

        resolved_username = _resolve_username(username)
        _execute_run(
            RunParams(
                env=env,
                pipeline_script=to_location(pipeline_script),
                output_dir=out_dir,
                launcher_cwd=launcher_cwd,
                username=resolved_username,
                analysis_dir=analysis_dir_val,
                resources_base=(
                    to_location(resources_base) if resources_base is not None else None
                ),
                upload_dst=upload_dst,
                file_lock=file_lock,
                rclone_options=rclone_options,
                name=name,
                singularity_cache=singularity_cache,
                callback=callback,
                work_dir_delete=work_dir_delete,
                project=project,
                disable=disable,
                enable=enable,
                analysis_name=analysis_name,
                driver_mode=driver_mode,
                passthrough_args=passthrough,
                run_sheet_content=rs_result.run_sheet_content,
                local_samplesheet=rs_result.local_samplesheet,
                env_overrides=env_override or [],
                pipeline=_effective_pipeline,
                pipeline_explicit=pipeline is not None and pipeline != "generic",
                wait=wait,
                wait_timeout_seconds=wait_timeout,
            )
        )
    except KeyboardInterrupt:
        logging.info("Cancelled by user")
        raise SystemExit(1)


@app.command()
def rclone(
    action: str = typer.Argument(help="rclone action (e.g. 'copy')."),
    src: str = typer.Argument(help="Source path or remote."),
    dst: str = typer.Argument(help="Destination path or remote."),
    num_partitions: int = typer.Option(
        10,
        help="Number of parallel rclone partitions.",
    ),
    work_dir: Optional[Path] = typer.Option(
        None,
        help="Directory to store partitions and rclone logs.",
    ),
    log_level: str = typer.Option(
        "INFO",
        help="Logging level.",
    ),
    include: Optional[list[str]] = typer.Option(
        None,
        help="Patterns to include for partitioning files.",
    ),
) -> None:
    """Parallel rclone copy using file partitioning."""
    level = getattr(logging, log_level.upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    if action != "copy":
        logging.error(
            "Unsupported rclone action: %s. Only 'copy' is supported.", action
        )
        raise SystemExit(1)

    from pipeline_launcher.rclone import rclone_copy

    rclone_copy(
        src=src,
        dst=dst,
        num_partitions=num_partitions,
        work_dir=work_dir,
        includes=include,
    )


if __name__ == "__main__":
    app()
