# Docker Guide

## Overview

Docker is a containerization platform that packages software applications and their dependencies into portable, lightweight containers.
This approach eliminates "it works on my machine" problems and ensures consistent execution across different computing environments, making it particularly valuable for bioinformatics workflows where reproducibility is critical.

This guide covers the basics of using Docker containers for bioinformatics analysis.
If you are already very comfortable using Docker for data analysis tasks, then you may not need this guide.

## Installing Docker

### Docker Installation

Install Docker on your system following the official documentation:

- **Linux**: Follow instructions at <https://docs.docker.com/engine/install/>
- **macOS**: Download Docker Desktop from <https://docs.docker.com/desktop/mac/>
- **Windows**: Download Docker Desktop from <https://docs.docker.com/desktop/windows/>

### Alternative Container Runtimes

XOOS Docker images are compatible with alternative container runtimes:

#### Podman

Podman can be used interchangeably with Docker commands:

```bash
# Replace 'docker' with 'podman' in any command
podman run --rm -u $(id -u):$(id -g) -v /path/to/data:/data -w /data image:tag command
```

#### Singularity/Apptainer

For HPC environments that use Singularity (now Apptainer), convert Docker images:

```bash
# Pull and convert Docker image to Singularity format
singularity pull docker://ghcr.io/roche-axelios/xoos/tool-name:latest

# Run with Singularity (note syntax differences)
singularity exec --bind /path/to/data:/data --pwd /data tool-name_latest.sif command

# Or use Singularity run
singularity run --bind /path/to/data:/data --pwd /data tool-name_latest.sif command
```

**Key Singularity differences:**

- Use `--bind` instead of `-v` for volume mounting
- Use `--pwd` instead of `-w` for working directory
- User ID mapping is handled automatically
- No need for `--rm` as containers don't persist by default

## Docker Best Practices for Bioinformatics

### Essential Docker Parameters

When running Docker containers for bioinformatics analysis, several parameters are essential:

```bash
docker run --rm -u $(id -u):$(id -g) -v /path/to/data:/data -w /data image:tag command
```

**Parameter explanations:**

- `--rm`: Automatically removes the container when it exits, preventing accumulation of stopped containers and saving disk space
- `-u $(id -u):$(id -g)`: Sets the user ID and group ID inside the container to match your host user, ensuring output files are owned by you rather than root
- `-v /host/path:/container/path`: Mounts host directories into the container, allowing access to input data and output results
- `-w /working/directory`: Sets the working directory inside the container where commands will be executed

### Example Usage

```bash
# Basic samtools analysis example
docker run --rm \
  -u $(id -u):$(id -g) \
  -v /home/user/data:/data \
  -v /home/user/results:/results \
  -w /data \
  staphb/samtools:1.17 \
  samtools view -b -q 20 /data/input.sam > /results/filtered.bam

# Multiple volume mounts for samtools workflow
docker run --rm \
  -u $(id -u):$(id -g) \
  -v /path/to/reference:/reference:ro \
  -v /path/to/input:/input:ro \
  -v /path/to/output:/output \
  -w /output \
  staphb/samtools:1.17 \
  samtools sort -@ 4 -o sorted.bam /input/unsorted.bam
```

## Accessing XOOS Docker Images

XOOS Docker images are hosted on GitHub Container Registry and can be accessed as follows:

### Pulling Images

```bash
# Pull the latest version
docker pull ghcr.io/roche-axelios/xoos/small_variant_caller:{{ space.vars.version }}
```

### Available Images

Browse available XOOS Docker images at:

- **GitHub Packages**: <https://github.com/orgs/Roche-AXELIOS/packages?repo_name=XOOS>
