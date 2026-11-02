# XOOS Pipeline Launcher

## Getting started

### Introduction

The XOOS Pipeline Launcher is a unified CLI tool (`xoos`) for launching Nextflow pipelines across different execution environments. It abstracts away the differences between local execution, Slurm HPC clusters, AWS Batch, and GCP Batch, providing a single interface for pipeline submission.

The launcher handles:

- Environment-specific configuration (container runtimes, scratch paths, job schedulers)
- Nextflow process management (JVM tuning, signal forwarding, log rotation)
- Artifact organization (output directories, work dir symlinks, problem reports)
- Post-run operations (result upload via rclone, work directory cleanup)

Pipeline-specific parameters (e.g. `--input`, `--outdir`, `--enable_subsample_run`) are passed through to Nextflow after a `--` separator, keeping the launcher CLI independent of any particular pipeline.

### Recommended system requirements

| Software          | Version  | Notes                                                                                                                                                                 |
|:------------------|:---------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Python            | >= 3.10  | Use operating system specific best practices.                                                                                                                         |
| Nextflow          | >= 25.10 | Follow steps outlined [in Nextflow documentation](https://www.nextflow.io/docs/stable/install.html#standalone-distribution).                                          |
| Java JRE          | 17-23    | Requirement of Nextflow. Use operating system specific best practices or see [Nextflow documentation](https://www.nextflow.io/docs/stable/install.html#requirements). |
| Container runtime | Latest   | Docker, Apptainer, Singularity, or Podman. See operating system and container runtime specific best practices.                                                        |

### Installation

```bash
pip install .
```

After installation the `xoos` command is available on PATH.

---

## Usage

### Running a pipeline

The `xoos run` command launches a Nextflow pipeline.
The `--env`, `--pipeline-script`, `--resources-base`, and `--analysis-dir` options are always required.
All pipeline-specific parameters are passed after the `--` separator.

#### AWS Batch

Run the xoos-nf-core germline pipeline on AWS Batch with a pre-existing run sheet.

Pass your own AWS Batch env config to `--env` (see [Environment Configuration](#environment-configuration)):

```bash
xoos run \
  --env /path/to/aws_batch.yaml \
  --pipeline-script /tmp/xoos-nf-core/main.nf \
  --resources-base s3://my-bucket/xoos-resources-1.1 \
  --analysis-dir s3://my-bucket/$USER/$(date +%Y%m%d)/my-analysis \
  -- \
  --input s3://my-bucket/run_sheet.csv \
  -profile germline_wgs
```

Run with auto-generated run sheet:

```bash
xoos run \
  --env /path/to/aws_batch.yaml \
  --pipeline-script /tmp/xoos-nf-core/main.nf \
  --resources-base s3://my-bucket/xoos-resources-1.1 \
  --analysis-dir s3://my-bucket/$USER/$(date +%Y%m%d)/my-analysis \
  --run-name small-test \
  --run-type SBX-D \
  --file-type basecall_fastq \
  --run-dir s3://my-bucket/test-data/small/fastq \
  --samplesheet s3://my-bucket/test-data/small/sample_sheet.csv \
  -- \
  -profile germline_wgs
```

The `--pipeline-script` can also be a cloud URI if the pipeline is already on S3:

```bash
xoos run \
  --env /path/to/aws_batch.yaml \
  --pipeline-script s3://my-bucket/xoos-nf-core-1.0/main.nf \
  --resources-base s3://my-bucket/xoos-resources-1.1 \
  --analysis-dir s3://my-bucket/my-analysis \
  -- \
  --input s3://my-bucket/run_sheet.csv \
  -profile germline_wgs
```

#### Slurm

Run on a Slurm HPC cluster. The `--project` parameter is required for Slurm environments to set the job comment.

Pass your own Slurm env config to `--env` (see [Environment Configuration](#environment-configuration)):

```bash
xoos run \
  --env /path/to/slurm.yaml \
  --pipeline-script /path/to/xoos-nf-core/main.nf \
  --resources-base /path/to/xoos-resources-1.1 \
  --analysis-dir /scratch/$USER/my-analysis \
  --project my-project \
  -- \
  --input /path/to/run_sheet.csv \
  -profile germline_wgs
```

Run with auto-generated run sheet:

```bash
xoos run \
  --env /path/to/slurm.yaml \
  --pipeline-script /path/to/xoos-nf-core/main.nf \
  --resources-base /path/to/xoos-resources-1.1 \
  --analysis-dir /scratch/$USER/my-analysis \
  --project my-project \
  --run-name run_001 \
  --run-type SBX-D \
  --file-type basecall_fastq \
  --run-dir /path/to/run_001 \
  --samplesheet /path/to/samples.csv \
  -- \
  -profile germline_wgs
```

#### Local

Run on a local machine with Docker:

```bash
xoos run \
  --env /path/to/local.yaml \
  --pipeline-script /path/to/xoos-nf-core/main.nf \
  --resources-base /path/to/xoos-resources-1.1 \
  --analysis-dir /path/to/my-analysis \
  -- \
  --input /path/to/run_sheet.csv \
  -profile germline_wgs
```

#### Notes

`--analysis-dir` is required and sets the root directory for all launcher artifacts (stage scripts, Nextflow state, work directory, and logs).
`--outdir` is a Nextflow pipeline parameter that is optional — it defaults to `{analysis-dir}/output` when not provided.

To add Nextflow profiles beyond those in the env config, pass `-profile` in the passthrough arguments — they are merged with the environment profiles automatically.

#### Run sheet generation

Instead of writing a run sheet CSV by hand and passing it via `--input`, the launcher can generate one from CLI flags.
When any of the five run sheet flags are provided, all five are required and `--input` must not appear in passthrough args.
For local and Slurm runs, the generated CSV is written to `<outdir>/run_sheet.csv` and `--input` is injected automatically.
For cloud runs, the CSV is uploaded to the stage directory in cloud storage.

| Flag            | Description                                                              | Valid values                                                   |
|:----------------|:-------------------------------------------------------------------------|:---------------------------------------------------------------|
| `--run-name`    | Unique run identifier. No spaces.                                        | String                                                         |
| `--run-type`    | Chemistry / instrument type.                                             | `SBX-D`, `SBX-FAST`, `YSU`                                    |
| `--file-type`   | Format and state of input data.                                          | `basecall_rdb`, `basecall_fastq`, `demultiplexed_fastq`, `bam`, `cram` |
| `--run-dir`     | Directory containing the run data.                                       | Directory path or cloud URI                                    |
| `--samplesheet` | Per-run samplesheet CSV (`sample_name`, `sample_sid`, optional columns). | File path or cloud URI                                         |

#### Execution environments

The `--env` flag selects an environment config JSON that determines the execution backend:

| Environment type | Description                                        |
|:-----------------|:---------------------------------------------------|
| local            | Runs Nextflow directly on the current machine      |
| slurm            | Generates a driver script and submits via `sbatch` |
| aws_batch        | Submits a containerized driver job to AWS Batch    |
| gcp_batch        | Submits a containerized driver job to GCP Batch    |

Pass the path to your environment config file via `--env`:

```bash
xoos run --env /path/to/my-cluster.yaml ...
```

See [Environment Configuration](#environment-configuration) for how to author one,
starting from the bundled `*.config.template` files.

#### Pipeline configuration

The `--pipeline` flag selects a pipeline-specific config that is composed with the env config before the run starts.
It injects the correct Nextflow profiles, config files, and pipeline params for the chosen pipeline without editing the env YAML.

The public xoosnf bundle omits site-specific overlays. Add an `env` mapping to a
custom pipeline config when your environment needs additional profiles or config
files.

| Value | Behavior |
|:------|:---------|
| *(omitted)* | Auto-detect: reads `manifest.name` from the pipeline script's sibling `nextflow.config`; resolves as a bundled name or falls back to generic |
| `xoosnf` | Load the bundled `xoosnf` pipeline config |
| `/path/to/custom.yaml` | Load a user-supplied pipeline config YAML |
| `generic` | Skip pipeline-specific config entirely; pass `--input` directly via passthrough args |

```bash
# Explicit bundled pipeline
xoos run --env /path/to/slurm.yaml --pipeline xoosnf ...

# External pipeline config
xoos run --env /path/to/slurm.yaml --pipeline /path/to/my-pipeline.yaml ...

# Generic — no pipeline config; supply --input yourself
xoos run --env /path/to/slurm.yaml --pipeline generic -- --input sheet.csv
```

The bundled pipeline config schema (`config.yaml`) accepts:

```yaml
pipeline_params:          # applied to all envs
  enable_low_gpu_memory_mode: true

env:                      # per-env overlays, keyed by env config stem
  slurm:
    profiles: [my_slurm_profile] # defined in the user-supplied config file
    config_files:
    - /path/to/my-slurm.config
```

#### Overriding environment config fields

The `--env-override` flag lets you override individual scalar fields in the loaded environment config without creating a separate YAML file. It is repeatable and uses dot-notation for nested fields.

```bash
xoos run \
  --env /path/to/slurm.yaml \
  --env-override driver.singularity_cache=/my/custom/cache \
  --env-override driver.scratch_base=/alt/scratch \
  ...
```

Overrides are applied after the YAML is loaded and parsed, so they take effect on the fully typed config object. Supported field types are `str`, `int`, `float`, `bool`, and their `Optional` variants. List and dict fields (e.g. `profiles`, `pipeline_params`) cannot be overridden — use a custom env YAML for those.

An empty value (`key=`) sets the field to an empty string. To clear an `Optional` field entirely, use the literal value `null`:

```bash
--env-override driver.singularity_cache=null
```

Common override examples:

| Override | Effect |
|:---------|:-------|
| `driver.singularity_cache=/my/path` | Use a custom Singularity image cache |
| `driver.scratch_base=/alt/scratch` | Change the scratch directory base |
| `aws_batch.region=eu-west-1` | Switch AWS region |
| `aws_batch.driver_vcpus=8` | Increase driver job CPU allocation |
| `gcp_batch.driver_machine_type=n2-standard-16` | Change GCP driver machine type |

### Parallel rclone copy

```bash
xoos rclone copy /local/data s3:bucket/path --num-partitions 10
```

Partitions files by size across N concurrent rclone processes for higher throughput.

---

## Environment Configuration

Environment configs are JSON files that describe how to execute a pipeline on a particular infrastructure. The `executor` field determines the backend.

### Config templates

The launcher ships starter Nextflow config templates under `nextflow_config/env/`,
one per execution backend. Copy the template matching your infrastructure,
replace the `<PLACEHOLDER>` values, and reference it from your environment config.

| Template                          | Executor  | Description                                  |
|:----------------------------------|:----------|:---------------------------------------------|
| `slurm.config.template`           | slurm     | Slurm HPC cluster                            |
| `standalone_server.config.template` | local   | Single server / local execution             |
| `aws_batch.config.template`       | aws_batch | AWS Batch                                    |
| `gcp_batch.config.template`       | gcp_batch | GCP Batch                                    |

Each backend also ships a matching environment config template under `env/`,
one per `*.config.template`. Copy it to `env/<name>.yaml`, replace the
`<PLACEHOLDER>` values (each field is documented inline), drop the `.template`
suffix, and select it with `--env <name>`.

| Env template                        | Executor  | Pairs with                          |
|:------------------------------------|:----------|:------------------------------------|
| `slurm.yaml.template`               | slurm     | `slurm.config.template`             |
| `standalone_server.yaml.template`   | local     | `standalone_server.config.template` |
| `aws_batch.yaml.template`           | aws_batch | `aws_batch.config.template`         |
| `gcp_batch.yaml.template`           | gcp_batch | `gcp_batch.config.template`         |

Driver script templates (Jinja2) live under `templates/`.

### Config fields

Scalar fields (string, int, float, bool) can be overridden from the CLI via `--env-override key=value` using dot-notation for nested fields. See [Overriding environment config fields](#overriding-environment-config-fields).

| Field             | Type     | Executor  | Description                                                                      |
|:------------------|:---------|:----------|:---------------------------------------------------------------------------------|
| `executor`        | string   | all       | `"local"`, `"slurm"`, `"aws_batch"`, or `"gcp_batch"`                            |
| `profiles`        | string[] | all       | Nextflow `-profile` values                                                       |
| `config_files`    | string[] | all       | Additional Nextflow config files (supports `@nextflow_config/` prefix)           |
| `pipeline_params` | object   | all       | Key-value pairs forwarded as `--key value` to Nextflow                           |
| `attributes`      | object   | slurm     | Slurm driver job settings; holds `driver_options` (see below)                    |
| `preamble`        | string   | slurm     | Shell lines prepended to the driver script (YAML block scalar; `string[]` accepted for legacy JSON) |
| `driver`          | object   | slurm     | `scratch_base` and `singularity_cache` paths; `cache_mode` (`shared`/`user`) sets the default cache strategy |
| `aws_batch`       | object   | aws_batch | AWS Batch submission parameters                                                  |
| `gcp_batch`       | object   | gcp_batch | GCP Batch submission parameters                                                  |

#### Slurm `attributes.driver_options`

`driver_options` is a free-form map of `#SBATCH` flags applied to the
**driver job** — the job that runs Nextflow itself. Each key becomes
`--<key>=<value>` in the generated `driver.sh`. Keys with `null`/empty
values are skipped.

| Key         | Example             | Driver `#SBATCH` directive    |
|:------------|:--------------------|:------------------------------|
| `partition` | `batch_cpu`         | `#SBATCH --partition=batch_cpu` |
| `qos`       | `<QOS>`             | `#SBATCH --qos=<QOS>`           |
| `account`   | `<ACCOUNT>`         | `#SBATCH --account=<ACCOUNT>`   |

> **Driver vs. task jobs.** `driver_options` only configures the driver
> job. The per-task sbatch submissions issued by Nextflow set their own
> `process.clusterOptions` (QOS, GRES) in the env's `.config`. By
> default, task jobs inherit the account from the driver job's
> environment. To charge task jobs to a specific account, set
> `slurm_account` in `pipeline_params` — it is forwarded as
> `--slurm_account` and consumed by `process.clusterOptions` via
> `params.getOrDefault('slurm_account', '')`. When unset, task jobs fall
> back to the inherited account.

Below is an example environment config for Slurm execution.

```json
{
  "executor": "slurm",
  "profiles": [
    "my_cluster"
  ],
  "attributes": {
    "driver_options": {
      "partition": "<CPU_PARTITION>",
      "qos": "<QOS>",
      "account": "<ACCOUNT>"
    }
  },
  "pipeline_params": {
    "slurm_account": "<ACCOUNT>"
  },
  "preamble": [
    "set +ue",
    "module load Java/17 Python/3",
    "set -ue",
    "",
    "export PATH=\"/path/to/nextflow:$PATH\""
  ],
  "config_files": [
    "/path/to/my_cluster.config"
  ],
  "driver": {
    "scratch_base": "/scratch/users/",
    "singularity_cache": "/path/to/singularity_cache"
  }
}
```

Copy `nextflow_config/env/slurm.config.template` to `my_cluster.config`,
replace the `<PLACEHOLDER>` values, and point `config_files` at it.

#### Bundled Nextflow configs

Config file paths prefixed with `@nextflow_config/` are resolved from the package's bundled `nextflow_config/` directory. This allows environment configs to reference hardware profiles shipped with the launcher without requiring absolute paths.

---

## Overview and CLI Options

### `xoos run` CLI options

Parameters in **bold** are required.

| Parameter                    | Description                                                                                                                      | Value(s)                                                                                |
|:-----------------------------|:---------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------|
| **--env**                    | Environment name or path to an env config JSON. Looked up as `env/{name}.json` in the package; falls back to a direct file path. | String or file path                                                                     |
| **--pipeline-script**        | Path or cloud URI to the Nextflow pipeline script (`.nf`). Cloud executors download the parent directory at runtime.             | File path or cloud URI (s3://, gs://)                                                   |
| **--resources-base**         | Path or cloud URI to the resources directory. Must contain a `resources.config` file.                                            | File path or cloud URI (s3://, gs://)                                                   |
| **--analysis-dir**           | Root directory or cloud URI for all pipeline artifacts (stage scripts, Nextflow state, work directory, and logs). `--outdir` defaults to `{analysis-dir}/output` when not provided. | File path or cloud URI (s3://, gs://)                             |
| --username                   | Override the username for cloud work directory paths and resource labels. Defaults to the current OS user.                       | String                                                                                  |
| --upload-dst                 | Destination for uploading results via rclone.                                                                                    | String (rclone remote path)                                                             |
| --file-lock                  | Path to a lock file for preventing concurrent executions.                                                                        | File path                                                                               |
| --rclone-options             | Options to pass to rclone (quoted string).                                                                                       | String [default: `""`]                                                                  |
| --name                       | Base name for analysis (used for job name and Nextflow `-name`).                                                                 | String                                                                                  |
| --singularity-cache          | Singularity cache strategy. Overrides the env config's `driver.cache_mode`.                                                       | `shared` or `user` [default: `driver.cache_mode`, else `shared`]                        |
| --callback                   | Command to invoke with START/COMPLETE/FAIL around the Nextflow run.                                                              | String (executable path)                                                                |
| --work-dir-delete            | Work directory cleanup policy.                                                                                                   | `delete-if-succeeded`, `delete-always`, `delete-never` [default: `delete-if-succeeded`] |
| --project                    | Project name for Slurm job comment. Required for HPC submission.                                                                 | String                                                                                  |
| --disable                    | Features to disable. Repeatable. Accepted values: `provenance-report`, `sample-report`, `per-sample-error-ignore`.              | String (repeatable)                                                                     |
| --enable                     | Features to enable. Repeatable. Accepted values: `pb-dynamic-scaling`.                                                           | String (repeatable)                                                                     |
| --env-override               | Override env config fields without editing the YAML. Repeatable. Uses dot-notation for nested fields (e.g. `driver.singularity_cache=/my/path`). See [Overriding environment config fields](#overriding-environment-config-fields). | `key=value` (repeatable) |

Any other arguments are forwarded to Nextflow as passthrough arguments.

#### Run sheet generation options

| Parameter     | Description                                                                   | Value(s)                                                       |
|:--------------|:------------------------------------------------------------------------------|:---------------------------------------------------------------|
| --run-name    | Run name for auto-generated run sheet. When set, all five flags are required. | String (no spaces)                                             |
| --run-type    | Run type.                                                                     | `SBX-D`, `SBX-FAST`, `YSU`                                    |
| --file-type   | Input file type.                                                              | `basecall_rdb`, `basecall_fastq`, `demultiplexed_fastq`, `bam`, `cram` |
| --run-dir     | Directory containing run data.                                                | Directory path or cloud URI                                    |
| --samplesheet | Per-run samplesheet CSV.                                                      | File path or cloud URI                                         |

#### Passthrough arguments

Any additional arguments after `--` are passed directly to Nextflow.
This allows users to specify any pipeline-specific parameters without the launcher needing to be aware of them.

The table below lists commonly used passthrough arguments, but any valid Nextflow argument can be used.

| Parameter | Description                                                                              | Value(s)                                |
|:----------|:-----------------------------------------------------------------------------------------|:----------------------------------------|
| --outdir  | Output directory for the pipeline. Defaults to `{analysis-dir}/output` when not provided. | Directory path or cloud URI             |
| --input   | Input run sheet.                                                                         | File path to a CSV with run information |
| -profile  | Nextflow profiles to activate. Can be specified multiple times.                          | String (profile name)                   |
| -c        | Additional Nextflow config file. Can be specified multiple times.                        | File path                               |

### `xoos rclone` CLI options

Parameters in **bold** are required.

| Parameter          | Description                                                    | Value(s)                          |
|:-------------------|:---------------------------------------------------------------|:----------------------------------|
| **action**         | rclone action.                                                 | `copy`                            |
| **src**            | Source path or remote.                                         | String                            |
| **dst**            | Destination path or remote.                                    | String                            |
| --num-partitions   | Number of parallel rclone workers.                             | Integer > 0 [default: `10`]       |
| --work-dir         | Directory for partition files and logs.                        | Directory path                    |
| --log-level        | Logging level.                                                 | String [default: `INFO`]          |
| --include          | File patterns to include.                                      | String (glob pattern)             |

---

## Appendix

### Output directory layout

The launcher organizes outputs beneath the `--analysis-dir` path:

```text
{analysis_dir}/
├── output/              Pipeline results (Nextflow --outdir)
├── nextflow/            Nextflow runtime artifacts (logs, trace, etc.)
│   └── work             Symlink to the actual Nextflow work directory
├── stage/               Batch driver scripts and staging artifacts
├── problem_reports/     Zipped logs collected on failure
└── .nextflow/           Nextflow state directory (used for -resume)
```

### Slurm driver mode

When submitting to Slurm, the launcher generates a driver script and submits it via `sbatch`. The driver script re-invokes `xoos run` with a `_driver` suffix on the environment name. The launcher detects this suffix and forces the local executor, so Nextflow runs directly on the allocated compute node instead of submitting another Slurm job.

### Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup, testing, and code style.
