# GCP Batch
<!-- markdownlint-disable MD024 -->

Run the XOOS pipeline on GCP Batch using the `gcp_batch` executor.

{% hint style="info" %}
For a full, step-by-step GCP deployment walkthrough (APIs, GCS bucket, service account,
networking, GPU quota, and a custom launcher environment), see
[Getting Started on GCP](../getting-started-gcp.md).
{% endhint %}

## Prerequisites

- A GCP project with billing enabled.
- `gcloud` CLI installed and authenticated.
- Batch API enabled on the project (`gcloud services enable batch.googleapis.com`).
- Docker images for the XOOS modules pushed to a container registry accessible from GCP (Artifact Registry or a public registry).

## Infrastructure

The following GCP resources are required to run the pipeline on GCP Batch:

1. **GCS bucket** for pipeline work directories, outputs, and staged configs.

2. **Service account** with the following roles:
   - `roles/batch.jobsEditor` — to create and manage Batch jobs.
   - `roles/storage.objectAdmin` — on the GCS bucket for read/write access.
   - `roles/iam.serviceAccountUser` — to allow the Batch job to run as this service account.

3. **Pipeline driver container image** — build and push a container image that includes Nextflow, Java, and the pipeline_launcher.
   This image runs the Nextflow driver process inside the GCP Batch job.

4. **Networking** — ensure the VPC and subnet used by Batch jobs have internet access (or access to your container registry and GCS).

## Environment config YAML

```yaml
executor: gcp_batch
profiles:
- my_gcp_profile
config_files:
- my_gcp_nextflow.config
gcp_batch:
  region: us-central1
  project_id: my-gcp-project
  gcs_bucket: my-pipeline-bucket
  pipeline_image_uri: us-docker.pkg.dev/my-project/my-repo/pipeline-driver:latest
  execution_service_account: pipeline-sa@my-gcp-project.iam.gserviceaccount.com
  cli_path: /usr/lib/google-cloud-sdk/bin/gcloud
  driver_machine_type: n2-standard-8
```

| Field                       | Description                                                          |
|:----------------------------|:---------------------------------------------------------------------|
| `region`                    | GCP region for Batch jobs and GCS.                                   |
| `project_id`                | GCP project ID.                                                      |
| `gcs_bucket`                | GCS bucket for work directories and staged files.                    |
| `pipeline_image_uri`        | Container image URI for the Nextflow driver job.                     |
| `execution_service_account` | Service account email used by Batch jobs.                            |
| `cli_path`                  | Path to `gcloud` inside the driver container.                        |
| `driver_machine_type`       | Machine type for the driver VM (e.g., `n2-standard-8`).              |

## Nextflow config profile

{% hint style="info" %}
Bundled GCP Nextflow config profiles are available (e.g., `gcp_batch_genomics_sbx`).
You can use a bundled profile via `config_files` in the env config, or supply your own via `-c` in passthrough args.
{% endhint %}

Your GCP Nextflow config should set `process.executor = 'google-batch'` and configure:

- Machine type selection based on resource requirements.
- Retry strategies for preemptible VM interruptions.
- Resource allocations per process label.

## Example run command

```bash
xoos run \
  --env /path/to/my_gcp_batch.yaml \
  --pipeline-script /path/to/xoosnf/main.nf \
  --resources-base gs://my-resources-bucket/xoos-resources-1.1 \
  --analysis-dir gs://my-pipeline-bucket/$USER/$(date +%Y%m%d)/my-analysis \
  -- \
  --input gs://my-pipeline-bucket/data/run_sheet.csv \
  -profile germline_wgs_duplex
```

## Detailed deployment guide

See [Getting Started on GCP](../getting-started-gcp.md) for a complete walkthrough of setting up
GCP resources and authoring a custom launcher environment.
