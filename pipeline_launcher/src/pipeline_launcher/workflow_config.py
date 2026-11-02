"""Parse workflow configuration YAML and derive pipeline parameters.

Translates a workflow configuration into the CLI arguments that
pipeline_launcher and xoos-nf-core expect.  The workflow config YAML is
an alternative to manually specifying --run-type, --file-type,
--samplesheet, etc.

Update the mapping tables when new analysis protocols are added.
"""

from __future__ import annotations

import csv
import logging
import re
from pathlib import Path

import yaml
from cloudpathlib import CloudPath
from pydantic import BaseModel, ConfigDict, Field, field_validator

from pipeline_launcher.location import LocationPath

# ---------------------------------------------------------------------------
# Pydantic models — enforce the WorkflowConfiguration schema
# ---------------------------------------------------------------------------

_SAMPLE_SID_RE = re.compile(r"^[ACGT]+$")
_SUPPORTED_SCHEMA_VERSIONS = {"0.4", "1.0"}


class Sample(BaseModel):
    """A multiplexed sample from the workflow configuration."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    sample_id: str = Field(alias="SampleId")
    index: str = Field(alias="Index")

    @field_validator("sample_id")
    @classmethod
    def _validate_sample_id(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("SampleId must not be empty")
        return v

    @field_validator("index")
    @classmethod
    def _validate_index(cls, v: str) -> str:
        if not _SAMPLE_SID_RE.fullmatch(v):
            raise ValueError(f"Index {v!r} must match [ACGT]+")
        return v


class GenomicResource(BaseModel):
    """A genomic resource file referenced by the analysis step."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    name: str = Field(alias="Name")
    type: str = Field(alias="Type")
    file_path: str = Field(alias="FilePath")
    sha256: str = Field(alias="Sha256")


class AnalysisConfig(BaseModel):
    """The Analysis section of a WorkflowConfiguration."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    device_name: str | None = Field(None, alias="DeviceName")
    device_software_version: str | None = Field(None, alias="DeviceSoftwareVersion")
    protocol_name: str | None = Field(None, alias="ProtocolName")
    protocol_configuration_name: str | None = Field(
        None, alias="ProtocolConfigurationName"
    )
    protocol_parameters: dict | None = Field(None, alias="ProtocolParameters")
    genomic_resource_files: list[GenomicResource] = Field(
        default_factory=list, alias="GenomicResourceFiles"
    )


class SequencingConfig(BaseModel):
    """The Sequencing section of a WorkflowConfiguration."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    device_name: str | None = Field(None, alias="DeviceName")
    device_software_version: str | None = Field(None, alias="DeviceSoftwareVersion")
    output_type: str | None = Field(None, alias="OutputType")
    protocol_name: str | None = Field(None, alias="ProtocolName")
    protocol_configuration_name: str | None = Field(
        None, alias="ProtocolConfigurationName"
    )
    protocol_parameters: dict | None = Field(None, alias="ProtocolParameters")


class SynthesisConfig(BaseModel):
    """The Synthesis section of a WorkflowConfiguration."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    device_name: str | None = Field(None, alias="DeviceName")
    device_software_version: str | None = Field(None, alias="DeviceSoftwareVersion")
    protocol_name: str | None = Field(None, alias="ProtocolName")
    protocol_configuration_name: str | None = Field(
        None, alias="ProtocolConfigurationName"
    )


class WorkflowConfig(BaseModel):
    """Parsed representation of a WorkflowConfiguration YAML."""

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    schema_version: str = Field(alias="SchemaVersion")
    workflow_order_name: str | None = Field(None, alias="WorkflowOrderName")
    samples: list[Sample] = Field(default_factory=list, alias="Samples")
    synthesis_tube_id: str | None = Field(None, alias="SynthesisTubeId")
    sequencing_tube_id: str | None = Field(None, alias="SequencingTubeId")
    xp_synthesis_completion_time: str | None = Field(
        None, alias="XpSynthesisCompletionTime"
    )
    sequencing: SequencingConfig | None = Field(None, alias="Sequencing")
    analysis: AnalysisConfig | None = Field(None, alias="Analysis")
    synthesis: SynthesisConfig | None = Field(None, alias="Synthesis")

    # Not part of the YAML — set after loading to track the source file.
    source_path: LocationPath | None = Field(None, exclude=True)


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------


def load_workflow_config(path: LocationPath) -> WorkflowConfig:
    """Load and validate a WorkflowConfiguration YAML.

    Raises ``ValueError`` on malformed YAML or schema validation errors.
    Logs a warning for unrecognised schema versions but does not fail.
    """
    try:
        raw = path.read_text()
    except FileNotFoundError:
        raise ValueError(f"Workflow config file not found: {path}") from None

    try:
        data = yaml.safe_load(raw)
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid YAML in {path}: {exc}") from exc

    if not isinstance(data, dict):
        raise ValueError(f"Expected YAML mapping in {path}, got {type(data).__name__}")

    try:
        config = WorkflowConfig.model_validate(data)
    except (ValueError, TypeError) as exc:
        raise ValueError(
            f"Workflow config validation failed for {path}: {exc}"
        ) from exc

    config.source_path = path

    # Warn on unrecognised schema versions but don't fail.
    # Normalise to major.minor so "0.4.0", "0.4.0-beta", and "0.4" all match.
    version_core = config.schema_version.split("-")[0]
    parts = version_core.split(".")
    major_minor = f"{parts[0]}.{parts[1]}" if len(parts) >= 2 else version_core
    if major_minor not in _SUPPORTED_SCHEMA_VERSIONS:
        logging.warning(
            "WorkflowConfiguration SchemaVersion %r is not in the "
            "supported set %s. Parsing will proceed but results may "
            "be incorrect.",
            config.schema_version,
            _SUPPORTED_SCHEMA_VERSIONS,
        )

    return config


# ---------------------------------------------------------------------------
# Derivation: file_type
# ---------------------------------------------------------------------------

# Alignment protocols — output is BAMs.
# Matched as prefixes: "BWA-based Alignment v2" matches "BWA-based Alignment".
# This accommodates version suffixes in the Everest spec protocol names.
_ALIGNMENT_PROTOCOLS = {"BWA-based Alignment", "Giraffe-based Alignment"}

# rdb2fastq converts RDB -> FASTQ; the basecall/consensus distinction
# is preserved from Sequencing.OutputType.
_FASTQ_CONVERSION_PROTOCOL = "rdb2fastq"

# Known keys in Analysis.ProtocolParameters that map to Nextflow params.
_KNOWN_PROTOCOL_PARAM_KEYS = {"alignmentOptions"}

# After rdb2fastq, RDB formats become their FASTQ equivalents.
# The sequencer only outputs RDB; FASTQ is produced by rdb2fastq.
_RDB2FASTQ_OUTPUT_TYPE = {
    "BASECALL_RDB": "basecall_fastq",
    "CONSENSUS_RDB": "demultiplexed_fastq",
}


def _is_alignment_protocol(protocol: str) -> bool:
    """Return True if *protocol* is a known alignment protocol."""
    return any(protocol.startswith(p) for p in _ALIGNMENT_PROTOCOLS)


def _file_type_from_rdb2fastq(ot: str | None) -> str:
    """Map Sequencing.OutputType to a FASTQ file_type after rdb2fastq."""
    if ot is None:
        raise ValueError(
            "Cannot derive file_type: rdb2fastq analysis ran but Sequencing.OutputType is absent"
        )
    ft = _RDB2FASTQ_OUTPUT_TYPE.get(ot)
    if ft is None:
        raise ValueError(
            f"Unexpected Sequencing.OutputType {ot!r} with rdb2fastq. "
            f"Expected one of: {sorted(_RDB2FASTQ_OUTPUT_TYPE)}"
        )
    return ft


def _file_type_from_output_type(ot: str | None) -> str:
    """Derive file_type from Sequencing.OutputType when no analysis ran."""
    if ot == "BASECALL_RDB":
        return "basecall_rdb"
    if ot == "CONSENSUS_RDB":
        raise ValueError(
            "Sequencing.OutputType is CONSENSUS_RDB but no rdb2fastq analysis was performed. "
            "CONSENSUS_RDB input is not supported without prior conversion to FASTQ."
        )
    if ot is None:
        raise ValueError(
            "Cannot derive file_type: neither Analysis.ProtocolName nor Sequencing.OutputType is present"
        )
    raise ValueError(
        f"Unexpected Sequencing.OutputType {ot!r}. The sequencer only outputs RDB formats; "
        f"FASTQ is produced by the rdb2fastq analysis protocol."
    )


def derive_file_type(config: WorkflowConfig) -> str:
    """Derive the pipeline ``file_type`` from the workflow configuration.

    Decision order:

    1. If the SAS ran an alignment protocol -> ``"bam"``
    2. If the SAS ran rdb2fastq -> map via ``Sequencing.OutputType``
    3. If no analysis section and ``BASECALL_RDB`` -> ``"basecall_rdb"``
    4. ``CONSENSUS_RDB`` without rdb2fastq -> error
    """
    protocol = config.analysis.protocol_name if config.analysis else None
    ot = config.sequencing.output_type if config.sequencing else None

    if protocol and _is_alignment_protocol(protocol):
        return "bam"

    if protocol and protocol.startswith(_FASTQ_CONVERSION_PROTOCOL):
        return _file_type_from_rdb2fastq(ot)

    # Warn if an unknown analysis protocol was present but didn't match.
    if protocol:
        logging.warning(
            "Unknown Analysis.ProtocolName %r; falling back to "
            "Sequencing.OutputType for file_type derivation.",
            protocol,
        )

    return _file_type_from_output_type(ot)


# ---------------------------------------------------------------------------
# Derivation: run_type
# ---------------------------------------------------------------------------

# Longest-prefix-first so more specific variants match before the fallback.
_PCN_PREFIX_TO_RUN_TYPE: list[tuple[str, str]] = [
    ("SBX-FAST", "SBX-FAST"),
    ("SBX-D", "SBX-D"),
]


def derive_run_type(config: WorkflowConfig) -> str:
    """Map ProtocolConfigurationName prefix to the pipeline's ``run_type``.

    Matches the longest known prefix of the PCN.  Falls back to ``"SBX-D"``
    with a WARNING when no prefix matches.
    """
    pcn = ""
    if config.sequencing:
        pcn = config.sequencing.protocol_configuration_name or ""
    for prefix, run_type in _PCN_PREFIX_TO_RUN_TYPE:
        if pcn.startswith(prefix):
            return run_type
    logging.warning(
        "No run_type mapping for ProtocolConfigurationName %r; defaulting to 'SBX-D'.",
        pcn,
    )
    return "SBX-D"


# ---------------------------------------------------------------------------
# Derivation: run_dir (parent of the YAML file)
# ---------------------------------------------------------------------------


def derive_run_dir(config: WorkflowConfig) -> LocationPath:
    """Return the parent directory of the workflow config YAML.

    Handles both local paths and cloud URIs.
    """
    if config.source_path is None:
        raise ValueError("Cannot derive run_dir: source_path not set")
    parent = config.source_path.parent
    if isinstance(parent, CloudPath):
        # CloudPath has no absolute(); return as-is since cloud URIs are
        # already absolute.
        return parent
    return Path(str(parent)).absolute()


# ---------------------------------------------------------------------------
# Samplesheet generation
# ---------------------------------------------------------------------------


def generate_samplesheet_csv(config: WorkflowConfig, out_dir: Path) -> Path:
    """Write a per-run samplesheet CSV from ``Samples[]``.

    Columns: ``sample_name``, ``sample_sid``
    Maps ``SampleId`` -> ``sample_name``, ``Index`` -> ``sample_sid``.

    Always writes to a local ``Path``.  For cloud executors, the caller
    is responsible for uploading the file (matching the existing
    ``_RunSheetResult.local_samplesheet`` pattern).

    Raises ``ValueError`` if no samples are present.
    """
    if not config.samples:
        raise ValueError("Cannot generate samplesheet: no samples in workflow config")

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "samplesheet.csv"
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_name", "sample_sid"])
        for s in config.samples:
            writer.writerow([s.sample_id, s.index])
    logging.info("Generated samplesheet from workflow config: %s", path)
    return path


# ---------------------------------------------------------------------------
# Derivation: pipeline params from Analysis.ProtocolParameters
# ---------------------------------------------------------------------------


def derive_pipeline_params(config: WorkflowConfig) -> dict[str, str]:
    """Map ``Analysis.ProtocolParameters`` to Nextflow ``--params``.

    Only emits params that differ from pipeline defaults.
    """
    params: dict[str, str] = {}
    if not config.analysis or not config.analysis.protocol_parameters:
        return params

    ap = config.analysis.protocol_parameters

    unknown = set(ap) - _KNOWN_PROTOCOL_PARAM_KEYS
    if unknown:
        logging.debug(
            "ProtocolParameters contains unmapped keys: %s. These will be ignored.",
            sorted(unknown),
        )

    # alignmentOptions.duplicateMarking
    dup = (ap.get("alignmentOptions") or {}).get("duplicateMarking")
    if dup == "Disable":
        params["dedup_strategy"] = "none"

    # alignmentOptions.penaltyForMismatch (BWA-specific)
    penalty = (ap.get("alignmentOptions") or {}).get("penaltyForMismatch")
    if penalty is not None:
        params["bwa_mismatch_penalty"] = str(penalty)

    return params


# ---------------------------------------------------------------------------
# Reference genome
# ---------------------------------------------------------------------------


def derive_reference_fasta(config: WorkflowConfig) -> str | None:
    """Extract the reference genome FASTA path from ``GenomicResourceFiles``."""
    if not config.analysis:
        return None
    for gr in config.analysis.genomic_resource_files:
        if gr.type == "REFERENCE_GENOME_FASTA":
            return gr.file_path
    return None


def warn_reference_integrity(
    config: WorkflowConfig, resources_base: LocationPath | None
) -> None:
    """Log a warning if the expected reference genome file is missing."""
    if not config.analysis:
        return
    for gr in config.analysis.genomic_resource_files:
        if gr.type != "REFERENCE_GENOME_FASTA":
            continue
        if resources_base is None:
            logging.warning(
                "Cannot verify reference genome %r: no --resources-base provided.",
                gr.name,
            )
            continue
        if isinstance(resources_base, CloudPath):
            logging.info(
                "Reference genome %r expected at %s/%s (cloud path, skipping existence check)",
                gr.name,
                resources_base,
                gr.file_path,
            )
            continue
        fasta_path = Path(str(resources_base)) / gr.file_path
        if not fasta_path.is_file():
            logging.warning(
                "Reference genome %r not found at %s",
                gr.name,
                fasta_path,
            )
