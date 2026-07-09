<!-- markdownlint-disable MD024 -->
# Workflow configuration YAML

The *workflow-configuration.yaml* file is written by the AXELIOS platform for each analysis task.
It describes the sequencing run, the multiplexed samples, and the analysis protocol that should be applied.
The pipeline_launcher reads this file and derives all run sheet parameters automatically, so you do not need to provide `--run-type`, `--file-type`, `--run-dir`, or `--samplesheet` separately.

## File location

The platform places *workflow-configuration.yaml* inside the run's output directory under `analyses/`:

```text
customer-storage/
└── analyses/
    └── <task-id>/
        └── workflow-configuration.yaml
```

Pass the full path to the pipeline_launcher with `--workflow-config`:

```bash
xoos run \
    --env /path/to/env.yaml \
    --pipeline-script /path/to/main.nf \
    --resources-base /path/to/resources \
    --analysis-dir /path/to/analysis \
    --workflow-config /customer-storage/analyses/task-001/workflow-configuration.yaml
```

## Top-level fields

| Field | Type | Required | Description |
|---|---|---|---|
| `SchemaVersion` | string | Yes | Version of the WorkflowConfiguration schema (e.g., `"0.4.0"`, `"1.0.0"`). |
| `WorkflowOrderName` | string | No | Human-readable name for the run order. Used as the run name in the pipeline. |
| `SynthesisTubeId` | string | No | Identifier for the synthesis tube. |
| `SequencingTubeId` | string | No | Identifier for the sequencing tube. |
| `Samples` | list | No | Array of multiplexed sample objects. See [Samples](#samples). |
| `Sequencing` | object | No | Sequencing protocol information. See [Sequencing](#sequencing). |
| `Analysis` | object | No | Analysis protocol and parameters. See [Analysis](#analysis). |
| `Synthesis` | object | No | Synthesis protocol information. See [Synthesis](#synthesis). |

### Schema version

The pipeline_launcher accepts schema versions `0.4` and `1.0`.
If the version is not recognized, processing continues with a warning.
Normalize the version to major.minor when checking compatibility (e.g., `"0.4.0"` and `"0.4.1"` both match `0.4`).

## Sequencing

The `Sequencing` section describes the output produced by the sequencing module.

| Field | Type | Required | Values | Description |
|---|---|---|---|---|
| `OutputType` | string | No | `BASECALL_RDB`, `CONSENSUS_RDB` | The format of the sequencer output files. |
| `ProtocolName` | string | No | — | Name of the sequencing protocol. |
| `ProtocolConfigurationName` | string | No | e.g., `"SBX-D v2"` | Name of the protocol configuration used during sequencing. |

**Output types:**

| `OutputType` | Description |
|---|---|
| `BASECALL_RDB` | Raw base-called reads in RDB format. Used for most germline workflows. |
| `CONSENSUS_RDB` | Duplex consensus reads in RDB format. Requires rdb2fastq conversion before alignment. |

```yaml
Sequencing:
  OutputType: BASECALL_RDB
  ProtocolName: SBX Sequencing Protocol v1
  ProtocolConfigurationName: SBX-D v2
```

## Samples

The `Samples` list contains one entry per multiplexed sample.
Each entry maps a sample name to the barcode index sequence used during library preparation.

| Field | Type | Required | Description |
|---|---|---|---|
| `SampleId` | string | Yes | Unique identifier for the sample. Must not be empty. |
| `Index` | string | Yes | Barcode index sequence. Must match the pattern `[ACGT]+` (uppercase bases only). |

The pipeline_launcher generates a samplesheet CSV from this list automatically.
The columns in the generated CSV are `sample_name` (from `SampleId`) and `sample_sid` (from `Index`).

```yaml
Samples:
  - SampleId: SampleA
    Index: ACGTACGT
  - SampleId: SampleB
    Index: TGCATGCA
  - SampleId: SampleC
    Index: AACCGGTT
```

## Analysis

The `Analysis` section specifies the analysis protocol and its configuration.

| Field | Type | Required | Description |
|---|---|---|---|
| `ProtocolName` | string | No | Name and version of the analysis protocol (e.g., `"BWA-based Alignment v1.1.0"`). |
| `ProtocolVersion` | string | No | Explicit version string for the protocol. |
| `ProtocolParameters` | object | No | Key-value parameters that override pipeline defaults. See [Protocol parameters](#protocol-parameters). |
| `GenomicResourceFiles` | list | No | Reference files used by the analysis. See [Genomic resource files](#genomic-resource-files). |

### Supported analysis protocols

| `ProtocolName` prefix | Derived `file_type` | Description |
|---|---|---|
| `BWA-based Alignment` | `bam` | Aligns reads with BWA. Produces a BAM file. |
| `Giraffe-based Alignment` | `bam` | Aligns reads with vg Giraffe. Produces a BAM file. |
| `rdb2fastq` | `basecall_fastq` or `demultiplexed_fastq` | Converts RDB to FASTQ. Output type depends on `Sequencing.OutputType`. |

The protocol name is matched as a prefix, so version suffixes such as `"BWA-based Alignment v2.0.1"` are recognized correctly.

### Protocol parameters

`ProtocolParameters` is a nested object.
Only the `alignmentOptions` key is mapped to pipeline parameters.

#### alignmentOptions

| Parameter | Values | Pipeline parameter | Description |
|---|---|---|---|
| `duplicateMarking` | `"Disable"`, `"Enable"` | `--dedup_strategy none` | When set to `"Disable"`, turns off duplicate marking. |
| `penaltyForMismatch` | integer | `--bwa_mismatch_penalty <n>` | Sets the BWA mismatch penalty score. |

```yaml
Analysis:
  ProtocolName: BWA-based Alignment v1.1.0
  ProtocolParameters:
    alignmentOptions:
      duplicateMarking: Disable
      penaltyForMismatch: 4
```

### Genomic resource files

Each entry in `GenomicResourceFiles` describes a reference file used during analysis.

| Field | Type | Required | Description |
|---|---|---|---|
| `Name` | string | Yes | Human-readable name for the resource (e.g., `"hg38"`). |
| `Type` | string | Yes | Resource type identifier. Use `REFERENCE_GENOME_FASTA` for reference genomes. |
| `FilePath` | string | Yes | Path to the file, relative to the resources base directory. |
| `Sha256` | string | Yes | SHA-256 checksum of the file for integrity verification. |

The pipeline_launcher extracts the path of the first `REFERENCE_GENOME_FASTA` entry and passes it as the reference genome to the pipeline.
If the resource file is not found at the expected path, the launcher logs a warning but continues.

```yaml
Analysis:
  GenomicResourceFiles:
    - Name: hg38
      Type: REFERENCE_GENOME_FASTA
      FilePath: genomes/hg38/hg38.fa
      Sha256: a3b4c5d6e7f8...
```

## Synthesis

The `Synthesis` section carries information about the synthesis step preceding sequencing.

| Field | Type | Required | Description |
|---|---|---|---|
| `ProtocolName` | string | No | Name of the synthesis protocol. |
| `ProtocolConfigurationName` | string | No | Name of the synthesis protocol configuration. |

```yaml
Synthesis:
  ProtocolName: SBX Synthesis Protocol v2
  ProtocolConfigurationName: SBX-S v2
```

## How pipeline parameters are derived

When you pass `--workflow-config`, the pipeline_launcher derives the following parameters from the YAML.
Explicit CLI flags always override derived values.

| Derived parameter | Source field(s) | Notes |
|---|---|---|
| `--run-name` | `WorkflowOrderName` | Falls back to the YAML filename stem if absent. |
| `--run-dir` | Path of the YAML file | The parent directory of *workflow-configuration.yaml*. |
| `--file-type` | `Analysis.ProtocolName` + `Sequencing.OutputType` | See the table below. |
| `--samplesheet` | `Samples[]` | A CSV is generated in the analysis output directory. |
| `--dedup_strategy` | `Analysis.ProtocolParameters.alignmentOptions.duplicateMarking` | Set to `none` when `Disable`. |
| `--bwa_mismatch_penalty` | `Analysis.ProtocolParameters.alignmentOptions.penaltyForMismatch` | Integer, forwarded as a string. |
| Reference FASTA | `Analysis.GenomicResourceFiles[type=REFERENCE_GENOME_FASTA].FilePath` | Resolved against `--resources-base`. |

### File type derivation

The `--file-type` parameter controls which input format the pipeline expects.
The launcher determines it in the following order:

| Condition | `file_type` |
|---|---|
| `Analysis.ProtocolName` starts with `"BWA-based Alignment"` or `"Giraffe-based Alignment"` | `bam` |
| `Analysis.ProtocolName` starts with `"rdb2fastq"` and `Sequencing.OutputType` is `BASECALL_RDB` | `basecall_fastq` |
| `Analysis.ProtocolName` starts with `"rdb2fastq"` and `Sequencing.OutputType` is `CONSENSUS_RDB` | `demultiplexed_fastq` |
| No `Analysis` section and `Sequencing.OutputType` is `BASECALL_RDB` | `basecall_rdb` |
| No `Analysis` section and `Sequencing.OutputType` is `CONSENSUS_RDB` | Error — conversion is required first. |

## Complete examples

### Minimal YAML (basecall RDB, no analysis)

A workflow that only specifies the sequencing output type and no downstream analysis.
The pipeline will receive raw basecall RDB files.

```yaml
SchemaVersion: "0.4.0"
Sequencing:
  OutputType: BASECALL_RDB
```

### BWA alignment workflow

A multiplexed germline workflow with three samples aligned to hg38.

```yaml
SchemaVersion: "0.4.0"
WorkflowOrderName: run-2024-03-15-batch-001
SynthesisTubeId: SYN-00123
SequencingTubeId: SEQ-00456
Samples:
  - SampleId: Sample001
    Index: ACGTACGT
  - SampleId: Sample002
    Index: TGCATGCA
  - SampleId: Sample003
    Index: AACCGGTT
Sequencing:
  OutputType: BASECALL_RDB
  ProtocolName: SBX Sequencing Protocol v1
  ProtocolConfigurationName: SBX-D v2
Analysis:
  ProtocolName: BWA-based Alignment v1.1.0
  ProtocolParameters:
    alignmentOptions:
      duplicateMarking: Enable
  GenomicResourceFiles:
    - Name: hg38
      Type: REFERENCE_GENOME_FASTA
      FilePath: genomes/hg38/hg38.fa
      Sha256: a3b4c5d6e7f8a1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6e7f8a1b2c3d4e5f6a7b8
Synthesis:
  ProtocolName: SBX Synthesis Protocol v2
  ProtocolConfigurationName: SBX-S v2
```

**Launch command:**

```bash
xoos run \
    --env /path/to/env.yaml \
    --pipeline-script /path/to/main.nf \
    --resources-base /customer-storage/resources \
    --analysis-dir /customer-storage/analyses/task-001 \
    --workflow-config /customer-storage/analyses/task-001/workflow-configuration.yaml
```

### FASTQ conversion from consensus RDB

A workflow that converts duplex consensus RDB output to demultiplexed FASTQ.

```yaml
SchemaVersion: "0.4.0"
WorkflowOrderName: duplex-run-2024-04-01
Samples:
  - SampleId: TumorSample
    Index: ACGT
  - SampleId: NormalSample
    Index: TGCA
Sequencing:
  OutputType: CONSENSUS_RDB
  ProtocolName: SBX Duplex Protocol v1
Analysis:
  ProtocolName: rdb2fastq v1.0.0
```

This configuration produces `demultiplexed_fastq` input for the pipeline.

### Alignment with duplicate marking disabled

Use this configuration when downstream tools perform their own duplicate handling.

```yaml
SchemaVersion: "1.0.0"
WorkflowOrderName: research-run-no-dedup
Samples:
  - SampleId: ResearchSample01
    Index: ACGTACGT
Sequencing:
  OutputType: BASECALL_RDB
  ProtocolConfigurationName: SBX-D v2
Analysis:
  ProtocolName: Giraffe-based Alignment v1.2.0
  ProtocolParameters:
    alignmentOptions:
      duplicateMarking: Disable
  GenomicResourceFiles:
    - Name: hg38-pangenome
      Type: REFERENCE_GENOME_FASTA
      FilePath: genomes/hg38-pangenome/hg38.fa
      Sha256: b1c2d3e4f5a6b7c8d9e0f1a2b3c4d5e6f7a8b9c0d1e2f3a4b5c6d7e8f9a0b1c2
```

## Overriding derived values

Explicit CLI flags take precedence over values derived from the workflow configuration YAML.
This allows you to correct or override individual fields without modifying the file.

```bash
xoos run \
    --env /path/to/env.yaml \
    --pipeline-script /path/to/main.nf \
    --resources-base /customer-storage/resources \
    --analysis-dir /customer-storage/analyses/task-001 \
    --workflow-config /customer-storage/analyses/task-001/workflow-configuration.yaml \
    --run-type SBX-D \
    --file-type bam
```

In this example, `--run-type` and `--file-type` override the values that would otherwise be derived from the YAML.

## Validation rules

The pipeline_launcher validates the YAML on load and exits with an error if validation fails.

| Rule | Error condition |
|---|---|
| `SchemaVersion` must be present | Missing field causes validation failure. |
| `Samples[].SampleId` must not be blank | Whitespace-only values are rejected. |
| `Samples[].Index` must match `[ACGT]+` | Any character outside `A`, `C`, `G`, `T` is rejected. |
| No extra top-level or nested fields are permitted | Unknown fields cause validation failure. |
| `CONSENSUS_RDB` without `rdb2fastq` analysis is not supported | The launcher returns an error. |
