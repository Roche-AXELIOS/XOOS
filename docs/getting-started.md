# Getting Started

## Philosophy

XOOS is a collection of secondary-analysis **modules** for SBX Duplex (SBX-D) and SBX Simplex (SBX-S) sequencing data.
The design is deliberately modular so you can adopt as much or as little as you need:

- **Run a module standalone.** Each module is a self-contained tool with its own CLI and Docker
  image. You can run any module on its own or wire it into an existing pipeline alongside your
  current tooling.
- **Run the end-to-end pipeline.** XOOS also ships a complete Nextflow pipeline (`xoosnf`) that
  chains the modules into a full secondary-analysis workflow you can run out of the box. The
  pipeline is written in an nf-core compliant manner, so it runs in any environment that can run
  Nextflow nf-core pipelines, and you can customize or extend it however you see fit. It can run
  on a standalone server, an HPC cluster, or in the cloud.

In short: use individual modules as building blocks, or use the pipeline as a ready-made,
extensible workflow.

## Modules

| Module | Description |
|:-------|:------------|
| [Demux](../demux/README.md) | Demultiplexes samples and trims adapters from raw reads; for duplex chemistries also performs consensus base calling. |
| [Read Collapser](../read_collapser/README.md) | Clusters reads by genomic position and/or UMI, then either marks duplicates or generates consensus reads. Supports germline WGS and target enrichment (TE). |
| [Pan-Genome Consensus Caller](../pangenome_consensus_caller/README.md) | Resolves duplex discordant bases using a pan-genome reference. |
| [Mitochondrial BAM Post-Process](../bam_post_process/README.md) | Recovers mitochondrial reads that align to nuclear mitochondrial DNA segments and assigns posterior mapping qualities. |
| [Alignment Metrics](../alignment_metrics/README.md) | Computes alignment and quality metrics from aligned SBX reads. |
| [Small Variant Caller](../small_variant_caller/README.md) | Filters and re-genotypes small variants (SNVs and short indels) using a machine-learning model. |
| [Copy Number Caller](../copy_number_caller/README.md) | Calls copy-number events and estimates purity/ploidy. |
| [STR Caller](../str_caller/README.md) | Genotypes and detects repeat expansions in short tandem repeats (STRs). |
| [Tumor Fraction Estimator](../tumor_fraction_estimator/README.md) | Estimates tumor fraction and detects sample contamination. |

## Pipeline Applications

The end-to-end pipeline supports the following applications, with chemistry-specific analysis profiles:

- **SBX-D Germline WGS**
- **SBX-S Germline WGS**
- **SBX-D Somatic Tumor/Normal WGS**
- **SBX-D Somatic Tumor-Only WGS**
- **SBX-D cfDNA WGS**

See [Applications](xoosnf/applications.md) for how to configure and run each one.

## Installation

### Installing a module

Installing a module is simple — pull its Docker image from the latest release:

```bash
docker pull ghcr.io/roche-axelios/xoos/<module>:latest
```

Browse available images and tags on the
[XOOS releases page](https://github.com/Roche-AXELIOS/XOOS/releases).

Each module image is self-contained and can be run standalone or dropped into an existing
pipeline. See the individual module pages under **Analysis Tools** for usage details.

### Installing the pipeline

The pipeline is distributed as a Python wheel. Download the wheel from the release assets and
install it into a virtual environment (replace `<version>` with the version you downloaded):

```bash
python -m venv .venv
source .venv/bin/activate
pip install pipeline_launcher-<version>-py3-none-any.whl
```

This installs the `xoos` launcher CLI. The pipeline is nf-core compliant and runs in any
environment that can run Nextflow nf-core pipelines.

The pipeline also needs a **resource bundle** (reference genome, indexes, BED files, and
trained models). Download it from the XOOS data portal on the Kamino website
([https://web.sbxdata.kamino.platform.navify.com/dashboard/](https://web.sbxdata.kamino.platform.navify.com/dashboard/))
and point the pipeline at it via `--resources-base`.

For the full setup — including downloading the resource bundle, creating an environment
configuration, and running your first analysis — see the pipeline's
[Getting started](xoosnf/README.md#getting-started) guide.
