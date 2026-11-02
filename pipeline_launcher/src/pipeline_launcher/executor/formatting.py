"""Shared command formatting and config generation for executor driver scripts.

Provides consistent quoted, backslash-continuation formatting across all
executor backends (AWS Batch, GCP Batch, Slurm), plus helpers for
generating Nextflow config snippets.

Values are single-quoted by default.  Values containing ``${`` are
double-quoted so the shell expands variable references (e.g.
``${PIPELINE_DIR}/main.nf`` when the driver sets ``PIPELINE_DIR`` at
runtime via ``mktemp -d``).
"""

from __future__ import annotations

import re
import shutil
import sys
from importlib import resources
from pathlib import Path
from typing import TYPE_CHECKING

from cloudpathlib import CloudPath
from jinja2 import Environment, FileSystemLoader, select_autoescape

from pipeline_launcher.location import to_location

if TYPE_CHECKING:
    from pipeline_launcher.executor.base import RunContext  # noqa: F401


def _single_quote(val: str) -> str:
    """Single-quote a value, escaping embedded single quotes."""
    return "'" + val.replace("'", "'\\''") + "'"


def _double_quote(val: str) -> str:
    """Double-quote a value, preserving shell variable expansion (``$var``, ``${var}``).

    Backslashes and double-quote characters in *val* are escaped so they
    are treated as literals; ``$`` is intentionally left unescaped so the
    shell expands variable references.
    """
    return '"' + val.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _quote(val: str) -> str:
    """Return *val* wrapped in the appropriate shell quotes.

    Values containing ``${`` are double-quoted so the shell expands the
    variable reference at runtime.  All other values are single-quoted.
    """
    if "${" in val:
        return _double_quote(val)
    return _single_quote(val)


def format_shell_command(args: list[str]) -> str:
    """Format an arbitrary shell command with backslash continuations.

    Same quoting rules as :func:`format_nf_command` but without
    prepending ``nextflow``.  The first element is the executable.
    Consecutive non-flag positional arguments (e.g. a subcommand such
    as ``run``) are placed on the same line as the executable, unquoted.
    """
    parts = _pair_args(args)
    if not parts:
        return ""

    # Collect the executable and any leading positional subcommand(s)
    # onto the first line, unquoted.  Stop at flags or at any positional
    # that looks like a path (contains "/" or "://"), which should be quoted.
    first_line: list[str] = []
    start = 0
    while start < len(parts) and not parts[start].startswith("-"):
        raw = parts[start].strip("'")
        if "/" in raw or "://" in raw:
            break
        first_line.append(raw)
        start += 1
    lines = [" ".join(first_line)]
    for part in parts[start:]:
        lines.append(f"    {part}")
    return " \\\n".join(lines)


def resolve_remote_path(
    local_path: Path | str | None,
    source_dir: Path | None,
    remote_base: str,
) -> str | None:
    """Map a local path to its remote equivalent inside the driver container.

    If *local_path* is relative to *source_dir*, returns the corresponding
    path under *remote_base*.  Otherwise returns the original path as a
    string (it may be a container-resident path).

    Returns ``None`` when *local_path* is ``None``.
    """
    if local_path is None:
        return None
    if source_dir is None:
        return str(local_path)
    try:
        rel = Path(local_path).relative_to(source_dir)
        return f"{remote_base}/{rel.as_posix()}"
    except ValueError:
        return str(local_path)


def build_nf_run_args(
    script: str,
    analysis_name_unique: str,
    *,
    config: str | None = None,
    extra_configs: list[str] | None = None,
    work_dir: str | None = None,
    output_dir: str | None = None,
    extra_args: list[str] | None = None,
) -> list[str]:
    """Build the Nextflow ``run`` argument list.

    Returns a list starting with ``"run"`` (no ``nextflow`` prefix).
    Callers should prepend ``"nextflow"`` as needed — either for subprocess
    use or before passing to :func:`format_shell_command` for formatted output.

    Assembles args in a consistent order: script, ``-name``, ``-resume``,
    user config, extra configs, ``-work-dir``, ``--outdir``, then any extra
    passthrough args.  ``-resume`` is always included.  ``--outdir`` is
    placed before ``extra_args`` so a passthrough ``--outdir`` overrides it
    (last-value-wins).
    """
    nf_run: list[str] = ["run", script, "-name", analysis_name_unique, "-resume"]
    if config is not None:
        nf_run.extend(["-c", config])
    for cfg_path in extra_configs or []:
        nf_run.extend(["-c", cfg_path])
    if work_dir is not None:
        nf_run.extend(["-work-dir", work_dir])
    if output_dir is not None:
        nf_run.extend(["--outdir", output_dir])
    if extra_args:
        nf_run.extend(a for a in extra_args if a != "-resume")
    return nf_run


def _sanitize_label(value: str) -> str:
    """Sanitize a value for use as a cloud resource label.

    GCP labels must be lowercase, max 63 chars, only ``[a-z0-9_-]``.
    AWS tags are more permissive but we normalize for consistency.
    """
    sanitized = re.sub(r"[^a-z0-9_-]", "-", value.lower())[:63]
    return sanitized.strip("-")


def _make_jinja_env() -> Environment:
    """Return a Jinja2 :class:`~jinja2.Environment` for the bundled templates.

    Auto-escaping is explicitly disabled via ``select_autoescape([])``
    because every template in this package renders a shell script, not
    HTML or XML.  Enabling HTML auto-escaping would corrupt shell syntax
    by escaping ``&``, ``<``, ``>``, and quotes.  The explicit call makes
    the intent clear and resolves the SonarQube
    "Disabling auto-escaping is security-sensitive" finding.
    """
    templates = Path(str(resources.files("pipeline_launcher") / "templates"))
    return Environment(
        loader=FileSystemLoader(str(templates)),
        autoescape=select_autoescape([]),
    )


def render_labels_config(user: str, analysis_name: str) -> str:
    """Render a Nextflow config snippet that sets ``process.resourceLabels``.
    Uses the deterministic analysis_name (no hash).
    """
    env = _make_jinja_env()
    template = env.get_template("resource_labels.config.jinja")
    return template.render(
        user=_sanitize_label(user),
        analysis_name=_sanitize_label(analysis_name),
    )


def _pair_args(args: list[str]) -> list[str]:
    """Group flags with their values and quote all values.

    Values are single-quoted unless they contain ``${``, in which case
    they are double-quoted to allow shell variable expansion.
    """
    parts: list[str] = []
    i = 0
    while i < len(args):
        arg = args[i]
        if (
            arg.startswith("-")
            and i + 1 < len(args)
            and not args[i + 1].startswith("-")
        ):
            # Flag with a value — keep on same line.
            parts.append(f"{arg} {_quote(args[i + 1])}")
            i += 2
        elif arg.startswith("-"):
            # Bare flag (e.g. -resume).
            parts.append(arg)
            i += 1
        else:
            # Positional value (e.g. script path).
            parts.append(_quote(arg))
            i += 1
    return parts


def extract_work_dir_from_args(args: list[str]) -> str | None:
    """Extract the value of ``-work-dir`` from a Nextflow argument list.

    Handles both the two-token form (``-work-dir VALUE``) and the
    equals form (``-work-dir=VALUE``).  Returns ``None`` when the flag
    is absent.
    """
    for i, arg in enumerate(args):
        if arg == "-work-dir" and i + 1 < len(args):
            return args[i + 1]
        if arg.startswith("-work-dir="):
            return arg.split("=", 1)[1]
    return None


# Flags consumed by the initial invocation that must not be forwarded to the
# driver script (run-sheet generation flags and --name are handled before the
# driver is created and must not trigger a second run-sheet generation).
_SKIP_FLAGS: frozenset[str] = frozenset(
    {
        "--name",
        "--run-name",
        "--run-type",
        "--file-type",
        "--run-dir",
        "--samplesheet",
    }
)

# Launcher flags whose values are local paths that must be canonicalized to
# absolute paths so the driver script works from a different cwd.
_CANONICALIZE_FLAGS: frozenset[str] = frozenset(
    {
        "--pipeline-script",
        "--resources-base",
    }
)

# Flags whose values are local file paths that must be resolved to absolute
# paths only when the value is an existing file (not a name like "xoosnf").
_CANONICALIZE_FILE_PATH_FLAGS: frozenset[str] = frozenset({"--env", "--pipeline"})


def _canonicalize_path_value(value: str) -> str:
    """Return *value* as-is if it is a cloud URI; resolve it to an absolute path otherwise."""
    loc = to_location(value)
    return value if isinstance(loc, CloudPath) else str(loc)


def _try_canonicalize_equals_form(arg: str) -> str:
    """Canonicalize the ``--flag=value`` form for flags in ``_CANONICALIZE_FILE_PATH_FLAGS``.

    Returns the canonicalized argument string when *arg* matches a flag in
    ``_CANONICALIZE_FILE_PATH_FLAGS`` using the equals form and the value is
    an existing local file path.  Returns ``arg`` unchanged in all other cases
    (flag matches but value is not a local file, or no flag matches).
    """
    for flag in _CANONICALIZE_FILE_PATH_FLAGS:
        if arg.startswith(f"{flag}="):
            value = arg.split("=", 1)[1]
            return f"{flag}={Path(value).resolve()}" if Path(value).is_file() else arg
    return arg


def _resolve_xoos_bin() -> str:
    """Return the absolute path to the xoos binary from sys.argv[0]."""
    argv0 = sys.argv[0]
    p = Path(argv0)
    if p.is_absolute():
        return str(p)
    resolved = shutil.which(argv0)
    return resolved if resolved else str(p.resolve())


def _collect_launcher_args() -> list[str]:
    """Extract launcher args from sys.argv, excluding passthrough and consumed flags.

    Reads sys.argv[2:] (skipping the executable and ``run`` subcommand).
    Stops consuming from sys.argv at ``--``; passthrough args are supplied
    separately by the caller.  Flags in ``_SKIP_FLAGS`` are dropped.
    Flags in ``_CANONICALIZE_FLAGS`` have their values canonicalized to
    absolute paths.  ``--env`` values that are file paths are resolved.
    """
    launcher_args: list[str] = []
    skip_next = False
    in_passthrough = False
    for i, arg in enumerate(sys.argv[2:], 2):
        if skip_next:
            skip_next = False
            continue
        if arg == "--":
            in_passthrough = True
            continue
        if in_passthrough:
            continue
        if arg in _SKIP_FLAGS:
            skip_next = True
            continue
        if arg in _CANONICALIZE_FLAGS:
            skip_next = True
            launcher_args.extend([arg, _canonicalize_path_value(sys.argv[i + 1])])
            continue
        if arg in _CANONICALIZE_FILE_PATH_FLAGS and Path(sys.argv[i + 1]).is_file():
            skip_next = True
            launcher_args.extend([arg, str(Path(sys.argv[i + 1]).resolve())])
            continue
        # Handle the --flag=value equals form; returns arg unchanged if no match.
        launcher_args.append(_try_canonicalize_equals_form(arg))
    return launcher_args


def build_driver_command(
    analysis_name: str,
    extra_args: list[str],
) -> list[str]:
    """Convert current sys.argv into a command for the driver script.

    Preserves the original --env value and adds --driver-mode so the
    re-invoked launcher uses the local executor.  The --name flag is
    replaced with --analysis-name to lock the run name.

    ``extra_args`` supplies the passthrough arguments (everything after
    ``--``) instead of extracting them from ``sys.argv``.  This is
    necessary because the initial invocation may inject ``--input``
    into the passthrough list after generating a run sheet — an
    argument that never appeared in the original ``sys.argv``.

    Relative paths in launcher-managed arguments are canonicalized to
    absolute paths since the driver re-invokes from a different cwd.
    """
    launcher_args = _collect_launcher_args()
    command = [_resolve_xoos_bin(), "run"]
    command.extend(launcher_args)
    command.extend(["--driver-mode", "--analysis-name", analysis_name, "--"])
    command.extend(extra_args)
    return command
