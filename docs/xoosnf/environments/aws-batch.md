# AWS Batch
<!-- markdownlint-disable MD024 -->

Run the XOOS pipeline on AWS Batch using the `aws_batch` executor.

{% hint style="info" %}
For a full, step-by-step AWS deployment walkthrough (Terraform/Terragrunt, SSM parameters, IAM,
and a custom launcher environment), see [Getting Started on AWS](../getting-started-aws.md).
{% endhint %}

## Prerequisites

- An AWS account with permissions to create and manage Batch, S3, and IAM resources.
- AWS CLI v2 installed and configured with a named profile or SSO.
- Docker images for the XOOS modules pushed to a container registry accessible from AWS Batch (ECR or a public registry).

## Infrastructure

The following AWS resources are required to run the pipeline on AWS Batch:

1. **S3 bucket** for pipeline work directories, outputs, and staged configs.
   The bucket must be in the same region as your Batch compute environments.

2. **IAM job role** with the following permissions:
   - `s3:GetObject`, `s3:PutObject`, `s3:ListBucket`, `s3:DeleteObject` on the S3 bucket.
   - `batch:DescribeJobs`, `batch:SubmitJob`, `batch:TerminateJob` for Nextflow to submit child jobs.
   - A trust policy allowing `ecs-tasks.amazonaws.com` to assume the role.

3. **Batch compute environments** — create one or more managed compute environments:
   - **CPU queues**: size them by instance type (e.g., `c7a.large` for small tasks, `c7a.8xlarge` for medium, `c7a.16xlarge` for high).
   - **GPU queue**: use GPU instances (e.g., `g6.16xlarge`) for alignment and variant calling processes that use Parabricks.
   - Set the compute environment to use Spot instances for cost savings, or On-Demand for guaranteed capacity.

4. **Batch job queues** — create a job queue for each compute environment, plus a **driver queue** for the Nextflow driver job (a small instance is sufficient, e.g., 2 vCPUs, 3 GB memory).

5. **Nextflow driver container image** — the launcher submits a containerized Nextflow driver job.
   The default image is `nextflow/nextflow:26.04.0`.
   Ensure this image is accessible from your Batch compute environment.

## Environment config YAML

Create a YAML file (or use a bundled one via `--env`):

```yaml
executor: aws_batch
profiles:
- my_aws_profile
pipeline_params: {}
config_files:
- my_aws_nextflow.config
aws_batch:
  region: us-west-2
  s3_bucket: my-pipeline-bucket
  driver_queue: my-driver-queue
  cli_path: /usr/local/bin/aws
  driver_image: nextflow/nextflow:26.04.0
  driver_vcpus: 2
  driver_memory_mib: 3072
  job_role_arn: arn:aws:iam::123456789012:role/my-nextflow-role
  expected_bucket_owner: '123456789012'
```

| Field                   | Description                                                                 |
|:------------------------|:----------------------------------------------------------------------------|
| `region`                | AWS region for Batch and S3.                                                |
| `s3_bucket`             | S3 bucket for work directories and staged files.                            |
| `driver_queue`          | Batch job queue for the Nextflow driver job.                                |
| `cli_path`              | Path to the AWS CLI inside the driver container.                            |
| `driver_image`          | Docker image for the Nextflow driver container.                             |
| `driver_vcpus`          | vCPUs allocated to the driver job.                                          |
| `driver_memory_mib`     | Memory (MiB) allocated to the driver job.                                   |
| `job_role_arn`          | IAM role ARN assumed by Batch jobs for S3 and Batch API access.             |
| `expected_bucket_owner` | AWS account ID that owns the S3 bucket (for ownership verification).        |

## Nextflow config profile

The environment config references a Nextflow config file via `config_files`.
This profile sets `process.executor = 'awsbatch'` and defines queue routing, retry strategies, and resource allocations.

{% hint style="info" %}
Config file paths prefixed with `@nextflow_config/` are resolved from the bundled `nextflow_config/` directory inside the pipeline_launcher package.
{% endhint %}

A typical AWS Batch Nextflow profile includes:

- `process.executor = 'awsbatch'` to route all tasks to Batch.
- A `queue` closure that selects the job queue based on CPU and memory requirements.
- Retry strategies for transient failures (exit codes 135-140, 255) and Spot interruptions.
- Resource allocations per process label (`process_single`, `process_low`, `process_medium`, `process_high`, `process_gpu`).

## Example run command

```bash
xoos run \
  --env /path/to/my_aws_batch.yaml \
  --pipeline-script /path/to/xoosnf/main.nf \
  --resources-base s3://my-resources-bucket/xoos-resources-1.1 \
  --analysis-dir s3://my-pipeline-bucket/$USER/$(date +%Y%m%d)/my-analysis \
  -- \
  --input s3://my-pipeline-bucket/data/run_sheet.csv \
  -profile germline_wgs_duplex
```

## Detailed deployment guide

See [Getting Started on AWS](../getting-started-aws.md) for a complete walkthrough of deploying
the multi-queue-batch Terraform module and authoring a custom launcher environment.
