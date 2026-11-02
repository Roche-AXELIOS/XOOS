"""Load and resolve environment configuration from YAML files.

Environment configs are YAML files that describe how to execute a Nextflow
pipeline on a particular infrastructure.  The ``executor`` field determines
which config type is constructed via Pydantic's discriminated union.

Resolution order for ``--env``:

  1. Bundled env/{name}.yaml inside the package
  2. env/{base}.yaml if name ends with "_driver" (strip suffix)
  3. Treat as a direct file path
  4. Fall back to env/default.yaml when --env is omitted

After loading, scalar fields can be overridden from the CLI via
``--env-override key=value`` (dot-notation for nested fields).
See :func:`apply_env_overrides`.
"""

from __future__ import annotations

import logging
import re
import types
from importlib import resources
from pathlib import Path
from typing import get_type_hints

import yaml
from cloudpathlib import CloudPath
from pydantic import BaseModel, TypeAdapter

from pipeline_launcher.config.model import EnvConfig, PipelineConfig, PipelineEnvOverlay

DRIVER_SUFFIX = "_driver"

_ENV_CONFIG_ADAPTER: TypeAdapter[EnvConfig] = TypeAdapter(EnvConfig)

_CONFIG_EXTENSION = ".yaml"


def is_driver_env(name: str) -> bool:
    return name.endswith(DRIVER_SUFFIX)


def strip_driver_suffix(name: str) -> str:
    return name[: -len(DRIVER_SUFFIX)] if is_driver_env(name) else name


def to_driver_env_name(base: str) -> str:
    return f"{base}{DRIVER_SUFFIX}"


def _bundled_env_dir() -> Path:
    """Return the path to the bundled env/ directory inside the package."""
    return Path(str(resources.files("pipeline_launcher") / "env"))


_BUNDLED_CONFIG_PREFIX = "@nextflow_config/"


def bundled_nfconfig(name: str) -> Path:
    """Return the path to a bundled Nextflow config file."""
    return Path(str(resources.files("pipeline_launcher") / "nextflow_config" / name))


def resolve_config_path(raw_path: str, launcher_cwd: Path | None = None) -> Path:
    """Resolve a config file path, handling the ``@nextflow_config/`` prefix.

    Paths prefixed with ``@nextflow_config/`` are resolved from the bundled
    ``nextflow_config/`` directory inside the pipeline_launcher package.
    Relative paths are resolved against *launcher_cwd* when provided.
    """
    if raw_path.startswith(_BUNDLED_CONFIG_PREFIX):
        return bundled_nfconfig(raw_path.removeprefix(_BUNDLED_CONFIG_PREFIX))
    path = Path(raw_path)
    if not path.is_absolute() and launcher_cwd is not None:
        path = launcher_cwd / path
    return path


def _find_config_file(env_dir: Path, stem: str) -> Path | None:
    """Return the config file for *stem* in *env_dir* if it exists."""
    candidate = env_dir / f"{stem}{_CONFIG_EXTENSION}"
    return candidate if candidate.exists() else None


def resolve_env_config_path(env: str | None, search_dir: Path | None = None) -> Path:
    """Resolve the environment config path from a name or file path.

    Looks for ``{name}.yaml`` in *search_dir* (or the bundled env
    directory).  *search_dir* overrides the bundled env directory,
    allowing callers to point at a project-specific env/ folder.
    """
    env_dir = search_dir if search_dir is not None else _bundled_env_dir()

    if env is None:
        default = _find_config_file(env_dir, "default")
        if default is None:
            logging.error(
                "No --env specified and no env/default.yaml found. "
                "Use --env to select an environment."
            )
            raise SystemExit(1)
        return default

    direct = _find_config_file(env_dir, env)
    if direct is not None:
        return direct

    if is_driver_env(env):
        base = _find_config_file(env_dir, strip_driver_suffix(env))
        if base is not None:
            return base

    as_path = Path(env)
    if as_path.exists():
        return as_path

    checked = [str(env_dir / f"{env}{_CONFIG_EXTENSION}")]
    if is_driver_env(env):
        base_stem = strip_driver_suffix(env)
        checked.append(str(env_dir / f"{base_stem}{_CONFIG_EXTENSION}"))
    checked.append(str(as_path))
    logging.error(
        f"Could not resolve env config for '{env}'. Checked: {', '.join(checked)}"
    )
    raise SystemExit(1)


def _load_raw(path: Path) -> dict:
    """Read a YAML file and return the parsed dict."""
    suffix = path.suffix.lower()
    if suffix != _CONFIG_EXTENSION:
        raise ValueError(
            f"Unsupported config file extension '{suffix}' for {path}. "
            f"Expected: {_CONFIG_EXTENSION}"
        )
    data = yaml.safe_load(path.read_text())
    # An empty file yields None from yaml.safe_load.
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(
            f"Expected a YAML mapping in {path}, got {type(data).__name__}"
        )
    return data


def load_env_config(path: Path) -> EnvConfig:
    """Parse a YAML env config file and return the typed config."""
    raw = _load_raw(path)
    # Default executor to "local" when omitted.
    raw.setdefault("executor", "local")
    return _ENV_CONFIG_ADAPTER.validate_python(raw)


# ---------------------------------------------------------------------------
# Pipeline config resolution and application
# ---------------------------------------------------------------------------

_MANIFEST_NAME_RE = re.compile(
    r'manifest\s*\{[^}]*name\s*=\s*[\'"]([^\'"]+)[\'"]',
    re.DOTALL,
)


def detect_pipeline_name(pipeline_source_dir: Path | CloudPath) -> str | None:
    """Extract manifest.name from a Nextflow config file via regex.

    Reads ``pipeline_source_dir / "nextflow.config"`` and extracts the
    pipeline name from the ``manifest`` block.  Accepts both local ``Path``
    and cloud ``CloudPath`` — both support ``/`` for joining and
    ``.read_text()`` for reading.  Returns ``None`` if the file is absent,
    unreadable, or the pattern is not found.  Never raises — all errors are
    logged at DEBUG level.

    Args:
        pipeline_source_dir: Directory containing the Nextflow pipeline source.

    Returns:
        The manifest name string, or ``None`` if it cannot be determined.
    """
    config_path = pipeline_source_dir / "nextflow.config"
    try:
        text = config_path.read_text()
    except Exception as exc:  # noqa: BLE001
        logging.debug("detect_pipeline_name: could not read %s: %s", config_path, exc)
        return None

    match = _MANIFEST_NAME_RE.search(text)
    if match is None:
        logging.debug("detect_pipeline_name: no manifest.name found in %s", config_path)
        return None

    return match.group(1)


def _bundled_pipeline_config(name: str) -> Path | None:
    """Return the bundled pipeline config path for *name* if it exists.

    Checks for ``pipeline/{name}/config.yaml`` inside the
    ``pipeline_launcher`` package's resources.  Returns ``None`` if the file
    does not exist.

    Args:
        name: Pipeline name to look up.

    Returns:
        Path to the bundled config YAML, or ``None`` if not found.
    """
    path = Path(
        str(resources.files("pipeline_launcher") / "pipeline" / name / "config.yaml")
    )
    return path if path.exists() else None


def resolve_pipeline_config_path(
    pipeline: str | None,
    pipeline_source_dir: Path | None,
    *,
    is_explicit: bool = True,
) -> Path | None:
    """Resolve the path to a pipeline-specific config YAML.

    Resolution order:

    1. ``pipeline == "generic"`` → return ``None`` (use generic behavior).
    2. ``pipeline is None`` → auto-detect from ``pipeline_source_dir``:
       a. If ``pipeline_source_dir`` is not ``None``, call
          :func:`detect_pipeline_name`.
       b. If a name is found, continue with that name to step 3.
       c. If not found, log INFO and return ``None`` (generic fallback).
    3. Try bundled: ``_bundled_pipeline_config(name)`` — return if found.
    4. Try as a direct file path: ``Path(pipeline)`` — return if it exists.
    5. Name not found: abort when ``is_explicit`` is ``True``; fall back to
       generic silently when ``False`` (auto-detected name not in bundle).

    Args:
        pipeline: Effective pipeline name/path, ``"generic"``, or ``None`` for
            auto-detection.
        pipeline_source_dir: Directory of the Nextflow pipeline source, used
            for auto-detection only.
        is_explicit: When ``True`` (default), an unresolvable name is a fatal
            error.  When ``False``, an unresolvable auto-detected name falls
            back to generic behavior silently.

    Returns:
        Resolved ``Path`` to the pipeline config YAML, or ``None`` when
        generic behavior should be used.

    Raises:
        SystemExit: When an explicit pipeline name cannot be resolved to a
            config file.
    """
    if pipeline == "generic":
        return None

    if pipeline is None:
        name: str | None = None
        if pipeline_source_dir is not None:
            name = detect_pipeline_name(pipeline_source_dir)
        if name is None:
            logging.debug(
                "Pipeline auto-detection: no manifest.name found, using generic behavior"
            )
            return None
        # Fall through with the detected name.
    else:
        name = pipeline

    # Step 3: try bundled.
    bundled = _bundled_pipeline_config(name)
    if bundled is not None:
        return bundled

    # Step 4: try as a direct file path.
    as_path = Path(name)
    if as_path.exists():
        return as_path

    # Step 5: name not found anywhere.
    if not is_explicit:
        logging.debug(
            "Auto-detected pipeline '%s' has no bundled config; using generic behavior",
            name,
        )
        return None
    checked = [
        str(resources.files("pipeline_launcher") / "pipeline" / name / "config.yaml"),
        str(as_path),
    ]
    logging.error(
        "Could not resolve pipeline config for '%s'. Checked: %s",
        pipeline,
        ", ".join(checked),
    )
    raise SystemExit(1)


def load_pipeline_config(path: Path) -> PipelineConfig:
    """Parse a pipeline config YAML and return a PipelineConfig.

    Args:
        path: Path to the pipeline config YAML file.

    Returns:
        Validated :class:`~pipeline_launcher.config.model.PipelineConfig`.

    Raises:
        ValueError: If the YAML does not conform to the ``PipelineConfig``
            schema (propagated from Pydantic).
    """
    raw = _load_raw(path)
    return PipelineConfig.model_validate(raw)


def apply_pipeline_config(
    env_config: EnvConfig,
    pipeline_config: PipelineConfig,
    env_name: str,
) -> None:
    """Mutate env_config in-place with pipeline-specific overlays.

    The ``_driver`` suffix is stripped from *env_name* before looking up the
    overlay so that a ``slurm_driver`` env resolves to the ``slurm`` overlay
    key.  Composition order for ``pipeline_params``:

    * env-level base params (lowest priority)
    * global pipeline params
    * env-specific overlay params (highest priority)

    ``profiles`` and ``config_files`` from the overlay are appended.
    Duplicate profiles are skipped (first-seen wins, order preserved).

    Args:
        env_config: The loaded environment config to mutate.
        pipeline_config: The loaded pipeline config supplying overrides.
        env_name: The raw env name (e.g. ``"slurm_driver"``); the
            ``_driver`` suffix is stripped before overlay lookup.
    """
    stripped = strip_driver_suffix(env_name)
    overlay = pipeline_config.env.get(stripped, PipelineEnvOverlay())

    env_config.pipeline_params = {
        **env_config.pipeline_params,
        **pipeline_config.pipeline_params,
        **overlay.pipeline_params,
    }

    # Append overlay profiles, skipping duplicates while preserving order.
    existing_profiles = set(env_config.profiles)
    for profile in overlay.profiles:
        if profile not in existing_profiles:
            env_config.profiles.append(profile)
            existing_profiles.add(profile)

    env_config.config_files.extend(overlay.config_files)


# ---------------------------------------------------------------------------
# Generic env-config overrides from CLI
# ---------------------------------------------------------------------------


def _coerce(value_str: str, target_type: type) -> str | int | float | bool | None:
    """Coerce a string value to *target_type*.

    Supports str, int, float, bool, and ``Optional[T]`` (where T is one of
    the above).  For any ``Optional[T]`` field, the literal string ``"null"``
    sets the value to ``None``.
    """
    # Unwrap Optional (Union[X, None]) to get the inner type.
    if isinstance(target_type, types.UnionType):
        args = target_type.__args__
        non_none = [a for a in args if a is not type(None)]
        if len(non_none) == 1:
            if value_str == "null":
                return None
            return _coerce(value_str, non_none[0])
        raise ValueError(
            f"Cannot override multi-type union {target_type}. "
            f"Only Optional[T] (single type | None) is supported."
        )

    if target_type is bool:
        if value_str.lower() in ("true", "1", "yes"):
            return True
        if value_str.lower() in ("false", "0", "no"):
            return False
        raise ValueError(f"Cannot convert {value_str!r} to bool")
    if target_type is int:
        return int(value_str)
    if target_type is float:
        return float(value_str)
    if target_type is str:
        return value_str

    raise ValueError(f"Unsupported target type {target_type} for override value")


def _resolve_target(config: BaseModel, key_path: str) -> tuple[BaseModel, str]:
    """Walk dot-notation *key_path* and return ``(parent_obj, field_name)``.

    Raises ``ValueError`` for unknown fields or non-model intermediates.
    """
    parts = key_path.split(".")

    obj: BaseModel = config
    for i, part in enumerate(parts[:-1]):
        if not isinstance(obj, BaseModel):
            raise ValueError(
                f"Cannot traverse into {'.'.join(parts[: i + 1])!r}: "
                f"{type(obj).__name__} is not a config model."
            )
        if not hasattr(obj, part):
            raise ValueError(
                f"Unknown field {'.'.join(parts[: i + 1])!r} on {type(obj).__name__}."
            )
        obj = getattr(obj, part)

    field_name = parts[-1]
    if not isinstance(obj, BaseModel):
        raise ValueError(
            f"Cannot set field on {'.'.join(parts[:-1])!r}: "
            f"{type(obj).__name__} is not a config model."
        )
    if not hasattr(obj, field_name):
        raise ValueError(f"Unknown field {key_path!r} on {type(obj).__name__}.")

    return obj, field_name


def _apply_single_override(config: BaseModel, key_path: str, value_str: str) -> None:
    """Resolve, coerce, and set a single override on *config*."""
    obj, field_name = _resolve_target(config, key_path)

    hints = get_type_hints(type(obj))
    if field_name not in hints:
        raise ValueError(
            f"Cannot determine type for {key_path!r} on {type(obj).__name__}."
        )

    target_type = hints[field_name]

    # Reject list/dict fields — only scalars are supported.
    origin = getattr(target_type, "__origin__", None)
    if origin in (list, dict):
        raise ValueError(
            f"Cannot override list/dict field {key_path!r}. "
            f"Only scalar fields are supported."
        )

    coerced = _coerce(value_str, target_type)
    setattr(obj, field_name, coerced)
    logging.info("env-override: %s = %r", key_path, coerced)


def apply_env_overrides(config: EnvConfig, overrides: list[tuple[str, str]]) -> None:
    """Apply pre-parsed ``(key, value)`` overrides to a loaded env config.

    Keys use dot-notation to address nested Pydantic model fields, e.g.
    ``driver.singularity_cache`` with value ``/my/path``.

    Only scalar fields (str, int, float, bool, and their Optional
    variants) can be overridden.  Attempting to set a list, dict, or
    unknown field raises ``ValueError``.

    Key format validation (identifier segments separated by dots) is
    handled by the CLI layer (``KeyValueParam``).
    """
    for key_path, value_str in overrides:
        _apply_single_override(config, key_path, value_str)
