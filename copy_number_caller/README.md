# Copy Number Caller

## Guide contents

- [Overview](#overview): module purpose, example applications, and system requirements.
- [Usage](docs/usage.md): command entry points and links to the full command guides.
  - [Call germline copy numbers](docs/usage-germline-wgs.md): germline copy number calling, inputs, outputs, tips, and CLI options.
  - [Call tumor-normal somatic copy numbers](docs/usage-tumor-normal-wgs.md): matched tumor-normal somatic calling, inputs, outputs, tips, and CLI options.
  - [Call tumor-only somatic copy numbers](docs/usage-tumor-only-wgs.md): tumor-only somatic calling, inputs, outputs, tips, and CLI options.
  - [Re-call somatic copy numbers](docs/usage-predict-somatic-cnv.md): re-calling copy numbers under a chosen purity/ploidy solution.
- [File formats](docs/file-formats.md): custom format schemas produced by Copy Number Caller.
- [Metrics](docs/metrics.md): summary metrics reported by each command.
- [Troubleshooting](docs/troubleshooting.md): common issues and resolutions.
- [Appendix](docs/appendix.md): shared concepts, terminology, allosome handling, and algorithm details.
- [Release notes](CHANGELOG.md): notable changes by version.

## Overview

The Copy Number Caller identifies genomic copy number variants by analyzing coverage patterns in DNA sequencing alignments.
Designed as a post-alignment processing tool, it enables the study of large-scale structural variants, from single-gene deletions to whole-chromosome aneuploidies.
It supports both germline copy number variant (GCNV) detection and allele-specific somatic copy number variant (SCNV) calling.

At a glance: **BAM + reference (+ tumor-normal VCF for matched somatic calling, or tumor VCF + population SNP VCF for tumor-only calling)** -> **segment (SEG) and VCF copy number callsets, plus log-ratio and B-allele fraction (BAF) visualization tracks**.

![High-level overview](docs/assets/high-level-overview.svg)

Copy Number Caller exposes four analysis commands.
Choose the command that matches your sample configuration; the [Usage](docs/usage.md) page walks through the decision in more detail.

- `germline-wgs` calls copy number variants from a germline normal (non-tumor) WGS sample.
- `tumor-normal-wgs` characterizes allele-specific somatic CNVs from a matched tumor-normal pair.
- `tumor-only-wgs` characterizes allele-specific somatic CNVs from a tumor sample when no matched normal is available.
- `predict-somatic-cnv` re-calls copy numbers from existing `tumor-normal-wgs` or `tumor-only-wgs` outputs under a chosen purity/ploidy solution, without re-running the full pipeline.

### Reference builds

hg38 is the only supported reference assembly. The bundled seed-segment, mappability, and blocklist resources are hg38-specific. If you attempt to use another assembly, you must provide your own matching resource files and ensure that the BAMs, VCF, and reference FASTA all use that assembly. See the [workflow guides](docs/usage.md) for command-specific requirements.

### Simplex and duplex base handling

The Copy Number Caller detects whether the input library is simplex or duplex to select library-specific analysis defaults, including the bundled mappability BigWig.
It derives its copy number signal (log-ratios) from per-interval read counts, and separately computes a per-base coverage depth that is used for GC correction and minimum-coverage filters.
Neither computation inspects base type, so simplex bases, including those in the simplex tail of a duplex read, are treated the same as duplex bases.

Library detection can inspect read-group metadata, the `YC` auxiliary tag, and base-quality patterns, but the copy number algorithms do not otherwise use base-type annotations.
The only read-level filtering applied is based on the SAM FLAG field through the `--exclude-flags` option (see the read filtering options for each subcommand).

### Example applications

Expand an example below to see the typical use case, the matching command guide, and a representative command.

<details>
<summary>Germline copy number variant detection</summary>

- **Data type:** Germline normal (non-tumor) WGS BAM.
- **Objective:** Detect inherited copy number gains and losses, from focal deletions to whole-chromosome aneuploidies.
- **Usage page:** [Call germline copy numbers example command](docs/usage-germline-wgs.md#example-command).
- **Command:** `copy_number_caller germline-wgs`

</details>

<details>
<summary>Matched tumor-normal somatic copy number alteration calling</summary>

- **Data type:** Matched tumor and normal WGS BAMs plus a tumor-normal VCF of germline SNPs.
- **Objective:** Call allele-specific somatic CNVs and estimate tumor purity and ploidy.
- **Usage page:** [Call tumor-normal somatic copy numbers example command](docs/usage-tumor-normal-wgs.md#example-command).
- **Command:** `copy_number_caller tumor-normal-wgs`

</details>

<details>
<summary>Tumor-only somatic copy number alteration calling</summary>

- **Data type:** Tumor WGS BAM, a single-sample tumor VCF, and a population SNP panel (for example, gnomAD).
- **Objective:** Call allele-specific somatic CNVs when no matched normal is available.
- **Usage page:** [Call tumor-only somatic copy numbers example command](docs/usage-tumor-only-wgs.md#example-command).
- **Command:** `copy_number_caller tumor-only-wgs`

</details>

<details>
<summary>Re-calling copy numbers under a forced purity/ploidy solution</summary>

- **Data type:** Segment, log-ratio, and BAF outputs from a prior `tumor-normal-wgs` or `tumor-only-wgs` run.
- **Objective:** Regenerate the SEG/VCF callset under a user-specified purity and ploidy without re-running the full pipeline.
- **Usage page:** [Re-call somatic copy numbers example command](docs/usage-predict-somatic-cnv.md#example-command).
- **Command:** `copy_number_caller predict-somatic-cnv`

</details>

### System requirements

|        | Requirements                                                                                   |
|:--------|:------------------------------------------------------------------------------------------------|
| CPU    | Modern processor with at least 16 cores. A processor supporting `AVX2` or `AVX-512` will result in much better performance. |
| Memory | At least 8 GiB.                                   |
