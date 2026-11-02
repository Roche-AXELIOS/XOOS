"""xoosnf run sheet plugin.

Wraps ``_generate_run_sheet_if_needed`` and ``_apply_workflow_config`` from
``pipeline_launcher.cli`` so the xoosnf pipeline participates in the standard
plugin discovery mechanism.
"""

from __future__ import annotations

from pathlib import Path

from pipeline_launcher.cli import (
    RunSheetParams,
    _apply_workflow_config,
    _generate_run_sheet_if_needed,
    _RunSheetResult,
)


class XoosnfRunSheetPlugin:
    """Run sheet plugin for the xoos-nf-core pipeline.

    Delegates to ``_apply_workflow_config`` when a workflow config YAML is
    provided, then calls ``_generate_run_sheet_if_needed`` to produce the
    CSV run sheet.
    """

    def generate(
        self,
        *,
        params: RunSheetParams,
        launcher_cwd: Path,
        is_cloud: bool,
        passthrough: list[str],
    ) -> _RunSheetResult:
        """Generate a run sheet for xoos-nf-core.

        If ``params.workflow_config`` is set, derives run sheet fields from
        the YAML via ``_apply_workflow_config``, then calls
        ``_generate_run_sheet_if_needed`` with the updated params.

        If no workflow config is set, calls ``_generate_run_sheet_if_needed``
        directly with the supplied params.

        Args:
            params: Run sheet parameters; may carry a ``workflow_config`` path.
            launcher_cwd: Working directory of the launcher process.
            is_cloud: True when the executor is a cloud executor.
            passthrough: Nextflow passthrough args (may be mutated in place).

        Returns:
            A _RunSheetResult with updated passthrough, optional run sheet
            content (cloud), and optional local samplesheet path.
        """
        run_name = params.run_name
        run_type = params.run_type
        file_type = params.file_type
        run_dir = params.run_dir
        samplesheet = params.samplesheet

        if params.workflow_config is not None:
            wf_result = _apply_workflow_config(
                params.workflow_config,
                run_name=run_name,
                run_type=run_type,
                file_type=file_type,
                run_dir=run_dir,
                samplesheet=samplesheet,
                resources_base=params.resources_base,
                launcher_cwd=launcher_cwd,
                passthrough=passthrough,
            )
            run_name = wf_result.run_name
            run_type = wf_result.run_type
            file_type = wf_result.file_type
            run_dir = wf_result.run_dir
            samplesheet = wf_result.samplesheet

        updated_params = RunSheetParams(
            run_name=run_name,
            run_type=run_type,
            file_type=file_type,
            run_dir=run_dir,
            samplesheet=samplesheet,
            workflow_config=params.workflow_config,
        )

        return _generate_run_sheet_if_needed(
            passthrough,
            updated_params,
            is_cloud=is_cloud,
        )


plugin = XoosnfRunSheetPlugin()
