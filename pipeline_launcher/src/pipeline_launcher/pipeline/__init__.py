"""Pipeline plugin protocol for run sheet and samplesheet generation."""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    from pipeline_launcher.cli import RunSheetParams, _RunSheetResult


class RunSheetPlugin(Protocol):
    """Protocol that all pipeline run sheet plugins must implement."""

    def generate(
        self,
        *,
        params: "RunSheetParams",
        launcher_cwd: Path,
        is_cloud: bool,
        passthrough: list[str],
    ) -> "_RunSheetResult":
        """Generate a run sheet and return the result.

        Args:
            params: Run sheet parameters (run_name, run_type, etc.).
            launcher_cwd: Working directory of the launcher process.
            is_cloud: True when the executor is a cloud executor.
            passthrough: Nextflow passthrough args (may be mutated).

        Returns:
            A _RunSheetResult with updated passthrough, optional content,
            and optional local samplesheet path.
        """
        ...


def load_run_sheet_plugin(pipeline_name: str) -> "RunSheetPlugin | None":
    """Import and return the plugin for a bundled pipeline, or None.

    Attempts to import ``pipeline_launcher.pipeline.<pipeline_name>`` and
    return its ``plugin`` attribute.  Returns ``None`` when the module does
    not exist or has no ``plugin`` attribute.

    Args:
        pipeline_name: The name of the bundled pipeline (e.g. ``"xoosnf"``).

    Returns:
        The pipeline's ``RunSheetPlugin`` instance, or ``None``.
    """
    try:
        mod = importlib.import_module(f"pipeline_launcher.pipeline.{pipeline_name}")
        return getattr(mod, "plugin", None)
    except ModuleNotFoundError:
        return None
