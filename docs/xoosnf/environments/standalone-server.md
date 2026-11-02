# Standalone server
<!-- markdownlint-disable MD024 -->

Run the XOOS pipeline on a single server using the `local` executor.

## Prerequisites

- A server with Docker installed (or Singularity/Apptainer).
- GPU drivers installed if using GPU-accelerated processes (Parabricks alignment and variant calling).
- Java 17+ and Nextflow installed.
- Python 3.10+ with the pipeline_launcher installed.

## Environment config YAML

A standalone server uses the `local` executor (or omits the `executor` field, which defaults to `local`):

```yaml
profiles:
- my_server
- docker
- gpu
config_files:
- my_server_nextflow.config
```

The `docker` and `gpu` profiles are defined in the xoosnf pipeline itself.
The `docker` profile enables Docker as the container runtime.
The `gpu` profile adds `--gpus all` to `docker.runOptions` and enables GPU support for Singularity/Apptainer.

## Nextflow config profile

Your server-specific Nextflow config should define resource limits matching the hardware:

```groovy
profiles {
    my_server {
        process {
            resourceLimits = [
                memory: 2000.GB,
                cpus: 256,
                time: 30.d
            ]

            // Resource allocations per process label
            withLabel:process_high {
                cpus   = 48
                memory = 120.GB
            }

            withLabel:process_gpu {
                ext.use_gpu = true
                accelerator = 1
                cpus = 64
                memory = 240.GB
            }
        }
    }
}
```

## Example run command

```bash
xoos run \
  --env /path/to/my_server.yaml \
  --pipeline-script /path/to/xoosnf/main.nf \
  --resources-base /path/to/xoos-resources-1.1 \
  --analysis-dir /data/analysis/my-run \
  -- \
  --input /data/run_sheet.csv \
  -profile germline_wgs_duplex
```
