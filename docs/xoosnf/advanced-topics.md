# Advanced topics

This guide covers advanced topics related to the xoosnf pipeline.

## Version pinning for reproducible runs

Controlled or offline environments can prevent automatic updates by pinning the pipeline source, container images, Nextflow configuration, and local caches before a run.
For container setup details, see the [Docker guide](../readme/docker-guide.md), which is published in GitBook at `/tutorials/docker-guide`.

### Pin container images

Use a versioned container tag for every XOOS tool instead of a floating tag such as `latest`.
For the strongest freeze, record the image digest after pulling the image and retain a local archive.
Set `XOOS_VERSION` once and reuse it for the container, cache, and pipeline source.

```bash
export XOOS_VERSION=1.1.0

docker pull "ghcr.io/roche-axelios/xoos/str_caller:${XOOS_VERSION}"
docker image inspect "ghcr.io/roche-axelios/xoos/str_caller:${XOOS_VERSION}" \
  --format '{{join .RepoDigests "\n"}}'
docker save --output "str_caller_${XOOS_VERSION}.tar" \
  "ghcr.io/roche-axelios/xoos/str_caller:${XOOS_VERSION}"
```

Docker stores images in the local daemon of the compute node.
If compute nodes clear the daemon between runs, restore the archive with `docker load --input "str_caller_${XOOS_VERSION}.tar"` before starting the pipeline.

On HPC or shared infrastructure, Apptainer or Singularity provides a user-managed image cache that can be retained between runs.
Configure the cache directory in the Nextflow profile that runs the pipeline so that the runtime uses the local image cache.
Use the executor scope that matches the selected profile:

```groovy
singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = "/data/containers/xoosnf-${System.getenv('XOOS_VERSION')}"
}
```

For a profile using Apptainer, configure the equivalent `apptainer` scope instead:

```groovy
apptainer {
    enabled    = true
    autoMounts = true
    cacheDir   = "/data/containers/xoosnf-${System.getenv('XOOS_VERSION')}"
}
```

Before switching to offline mode, run the workflow once while online with the same profile and cache directory.
Nextflow then resolves each image and stores it using the filename it will look up later.
A manual `apptainer pull` in the current directory is not sufficient because its filename may not match Nextflow's URL-derived cache filename.

### Pin the pipeline source

Clone the `xoosnf` repository at the release tag that the run must use, and record the commit ID for the run record.

```bash
git clone --branch "${XOOS_VERSION}" --depth 1 \
  https://github.com/Roche-Axelios/xoosnf "xoosnf-${XOOS_VERSION}"
git -C "xoosnf-${XOOS_VERSION}" rev-parse HEAD
```

Run the pinned checkout directly with Nextflow or pass its local `main.nf` to the pipeline launcher.

```bash
nextflow run "/your/local_path/xoosnf-${XOOS_VERSION}/main.nf" ...

xoos run \
  --env my_environment \
  --pipeline-script "/your/local_path/xoosnf-${XOOS_VERSION}/main.nf" \
  ...
```

### Disable online lookups and container pulls

Set `NXF_OFFLINE` before running Nextflow to block remote dependency and GitHub lookups.

```bash
export NXF_OFFLINE=true
```

For Docker, add `--pull=never` to the profile's existing `docker.runOptions` when every required image is already loaded on the compute node.
This Docker option prevents the runtime from attempting a registry pull during an offline run.

```groovy
docker {
    enabled    = true
    runOptions = '--pull=never' // Append to any existing run options.
}
```

For Apptainer or Singularity, use the local `cacheDir` shown above and verify that every process resolves to a cached image before starting a run.

### Verify the versions used by a run

Inspect the **Software Versions** section of the final MultiQC report to confirm the tool versions that executed.
Keep the report together with the pipeline commit ID, container digests or image archive checksums, and the Nextflow configuration used for the run.

## Performance tuning

### Reading the execution report and timeline

Every pipeline run produces four report files under `{outdir}/pipeline_info/`:

| File | Purpose |
|:-----|:--------|
| `execution_report_*.html` | Per-task resource usage: CPU, memory, duration, and read/write I/O. |
| `execution_timeline_*.html` | Gantt chart showing when each task started, ran, and finished. |
| `execution_trace_*.txt` | Tab-separated trace data (machine-readable version of the report). |
| `pipeline_dag_*.html` | Directed acyclic graph of the pipeline's task dependencies. |

Open the **execution report** in a browser.
The key metric to look at is the **% allocated** column for CPU and memory.
This shows how much of the requested resource each task actually used.

- **CPU % allocated well below 100%** means the task was given more CPUs than it can use.
  Reduce `cpus` for that process label or process name in your Nextflow config.
- **Memory % allocated well below 100%** means the task was given more memory than its peak RSS.
  Reduce `memory` for that process.
  Leave at least 20% headroom above peak RSS to account for run-to-run variance.
- **Memory % allocated near or above 100%** means the task is at risk of being OOM-killed.
  Increase `memory` or add a retry multiplier (e.g., `memory = { 64.GB * task.attempt }`).
- **Duration** values that are much shorter than the `time` limit indicate the time limit can be tightened.

The **execution timeline** shows the overall pipeline wall time and parallelism.
Look for:

- **Long sequential gaps** — tasks that block the pipeline because nothing else can run in parallel.
  These are candidates for splitting into smaller chunks or increasing parallelism.
- **Tasks that dominate wall time** — focus optimization efforts on these.
- **Idle periods** — may indicate that the executor queue size is too small or that resource limits are preventing tasks from being scheduled.

### Updating your environment Nextflow config

After reviewing the execution report, update the resource allocations in your environment's Nextflow config profile.
The config uses process labels and process names to set resources:

```groovy
profiles {
    my_environment {
        process {
            // Default resources for all tasks
            cpus   = { 1 * task.attempt }
            memory = { 4.GB * task.attempt }
            time   = { 4.h * task.attempt }

            // Override by process label
            withLabel:process_high {
                cpus   = 32
                memory = { 64.GB * task.attempt }
                time   = { 12.h * task.attempt }
            }

            // Override by specific process name
            withName: PARABRICKS_GIRAFFE_SINGLE_END_MULTIPLE_FASTQS {
                cpus   = 48
                memory = 128.GB
                time   = 16.h
            }
        }
    }
}
```

The `* task.attempt` pattern increases the resource on each retry attempt, which helps tasks that fail due to transient resource pressure.

{% hint style="info" %}
Use the `execution_trace_*.txt` file for scripted analysis.
It is a TSV file with columns for `name`, `hash`, `status`, `exit`, `realtime`, `%cpu`, `peak_rss`, `peak_vmem`, `rchar`, `wchar`, and more.
You can import it into a spreadsheet or process it with `awk`/`pandas` to identify systematic over- or under-allocation.
{% endhint %}

---

## Custom process overrides

You can override any Nextflow process directive (command-line arguments, container image, memory, CPU, time) without modifying the pipeline source code.
Create a custom Nextflow config file and pass it to the pipeline_launcher via `-c` in passthrough args.

### Overriding resource allocations

```groovy
// my_overrides.config
process {
    withName: XOOS_SMALL_VARIANT_CALLER {
        cpus   = 16
        memory = 48.GB
        time   = 8.h
    }
}
```

```bash
xoos run \
  --env my_environment \
  --pipeline-script /path/to/xoosnf/main.nf \
  --resources-base /path/to/resources \
  --analysis-dir /path/to/analysis \
  -- \
  --input run_sheet.csv \
  -profile germline_wgs_duplex \
  -c my_overrides.config
```

### Overriding the container image

To use a different version of a module's container image:

```groovy
// my_overrides.config
process {
    withName: XOOS_SMALL_VARIANT_CALLER {
        container = 'my-registry.example.com/xoos/small_variant_caller:2.0.0'
    }
}
```

### Overriding command-line arguments

Pipeline processes use `ext.args` to pass additional command-line arguments to the underlying tool.
You can append or replace arguments:

```groovy
// my_overrides.config
process {
    withName: XOOS_ALIGNMENT_METRICS {
        ext.args = '--min-mapping-quality 10 --min-base-quality 20'
    }
}
```

{% hint style="warning" %}
Overriding `ext.args` replaces the entire argument string for that process.
Check the process's `.nf` file to see what arguments are set by default before overriding.
{% endhint %}

### Overriding error strategy

To change how a specific process handles failures:

```groovy
// my_overrides.config
process {
    withName: XOOS_DEMUX {
        errorStrategy = 'terminate'
        maxRetries    = 0
    }
}
```

### Multiple overrides in one file

You can combine multiple `withName` blocks in a single config file:

```groovy
// my_overrides.config
process {
    withName: XOOS_SMALL_VARIANT_CALLER {
        cpus   = 16
        memory = 48.GB
    }

    withName: XOOS_ALIGNMENT_METRICS {
        memory = 32.GB
        time   = 4.h
    }

    withName: PARABRICKS_GIRAFFE_SINGLE_END_MULTIPLE_FASTQS {
        container = 'my-registry.example.com/xoos/alignment:custom'
        memory    = 256.GB
    }
}
```

{% hint style="info" %}
Config files passed via `-c` are applied after all profile configs.
This means your overrides take precedence over the environment and pipeline profile settings.
{% endhint %}

---

## Resume and retries

### How `-resume` works

Nextflow caches the result of every task in the work directory.
The pipeline_launcher automatically passes `-resume` to Nextflow on every run, so re-running a pipeline with the same `--analysis-dir` will reuse cached results.
Nextflow checks whether each task's inputs, script, and container image match a previous execution.
If they match, the cached result is reused and the task is not re-executed.

```bash
# Re-run with the same --analysis-dir to resume from cached results
xoos run \
  --env my_environment \
  --pipeline-script /path/to/xoosnf/main.nf \
  --resources-base /path/to/resources \
  --analysis-dir /path/to/analysis \
  -- \
  --input run_sheet.csv \
  -profile germline_wgs_duplex
```

Resume requires:

- The **work directory** from the previous run must still exist.
  If `--work-dir-delete` is set to `delete-if-succeeded` (the default), the work directory is deleted after a successful run and resume is not possible.
  Use `--work-dir-delete delete-never` if you plan to resume.
- The **`.nextflow/` state directory** must be present at `{analysis-dir}/.nextflow/`.
  This directory contains the execution history and cache metadata.
- The **same `--analysis-dir`** must be used so Nextflow can find the previous state.

A task is re-executed when any of the following change:

- Input files (content or path).
- The process script (`.command.sh`).
- The container image.
- Process directives that affect execution (CPU, memory, etc. do **not** invalidate the cache — only the script and inputs do).

{% hint style="info" %}
Changing resource allocations (CPU, memory, time) does **not** invalidate the cache.
You can increase memory for a failing task and resume without re-running tasks that already succeeded.
{% endhint %}

### How task-level retries work

Nextflow retries tasks automatically based on the `errorStrategy` and `maxRetries` directives.
The xoosnf pipeline configures retries at two levels:

#### Base retry strategy (all environments)

Defined in `conf/base.config`:

```groovy
errorStrategy = { task.exitStatus in ((130..145) + 104 + 175) ? 'retry' : 'finish' }
maxRetries    = 1
```

This retries tasks that fail with exit codes in the range 130–145 (signals like SIGKILL, SIGTERM), 104 (connection reset), or 175 (memory limit).
All other failures cause the pipeline to finish remaining tasks and then exit.

#### Environment-specific retry strategies

Environment Nextflow configs override the base strategy with more retries and environment-specific exit code handling.
All environment configs set `maxRetries = 3`.

**Slurm configs:**

- Retry on exit codes 135–140, 125, 255, and `Integer.MAX_VALUE`.
- Exit code 125 (Docker/container startup failure) is retried because container pulls on shared filesystems can fail transiently.
- `Integer.MAX_VALUE` retries handle Slurm preemption, but without backoff — the task is resubmitted immediately.

**AWS Batch Spot config:**

- Retry on exit codes 135–140 and 255.
- `Integer.MAX_VALUE` (Spot interruption) is retried with exponential backoff: the task sleeps `2^attempt * 10` seconds before resubmitting.
- On Spot retries, the `spotRetry` flag is set (`task.ext.spotRetry = true`).
  This tells the Parabricks dynamic scaling config to preserve the previous task's resource allocation instead of scaling up, since the failure was caused by instance reclamation, not insufficient resources.

**AWS Batch On-Demand config:**

- Retry on exit codes 135–140 and 255.
- Does not handle `Integer.MAX_VALUE` because On-Demand instances are not subject to Spot interruptions.

| Config | Retryable exit codes | Spot/preemption handling | Backoff |
|:-------|:---------------------|:------------------------|:--------|
| `base.config` | 130–145, 104, 175 | None | None |
| Slurm | 135–140, 125, 255, `MAX_VALUE` | Immediate retry | None |
| AWS Batch Spot | 135–140, 255, `MAX_VALUE` | Retry + `spotRetry` flag | `2^attempt * 10s` |
| AWS Batch On-Demand | 135–140, 255 | None | None |

#### Resource scaling on retry

Many resource allocations use `task.attempt` to increase resources on each retry:

```groovy
memory = { 64.GB * task.attempt }
time   = { 12.h * task.attempt }
```

On the first attempt, memory is 64 GB.
On the second attempt (first retry), it increases to 128 GB.
On the third attempt, 192 GB.
This helps tasks that fail due to data-dependent memory spikes.

#### Per-sample error handling

By default, the pipeline uses `ignore` as the fallback error strategy for non-retryable failures.
This means a single sample's failure does not stop the pipeline — other samples continue processing.
The launcher injects `workflow.failOnIgnore = true` so the pipeline still exits with a non-zero code when any task was ignored.

To switch to strict mode where any failure stops the pipeline immediately:

```bash
xoos run ... --disable per-sample-error-ignore
```

This changes the fallback error strategy from `ignore` to `finish`.
