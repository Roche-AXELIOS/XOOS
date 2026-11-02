# process-mitochondrial-bam file formats

This page describes the text formats used by `process-mitochondrial-bam`.
Usage details are on the [Recruiting reads](usage-recruit.md) and [Assigning mitochondrial reads](usage-call.md) pages.

## NuMT anchor TSV

Applies to: `process-mitochondrial-bam recruit`, `process-mitochondrial-bam call`.

The tab-separated file contains one polymorphic NuMT anchor per row and has no header.

| Column | Name | Description |
|:--------|:------|:-------------|
| 1 | Chromosome | Nuclear reference contig. |
| 2 | Position | Zero-based nuclear insertion anchor. |
| 3 | Insertion length | NuMT insertion length in bp; optional for `recruit` and required for the `call` length model. |

Back to [Recruiting reads input](usage-recruit.md#input) and [Assigning mitochondrial reads input](usage-call.md#input).

## Recruited read-name list

Applies to: `process-mitochondrial-bam recruit`.

The file contains one read name per line.
Names are deduplicated and sorted.

Back to [Recruiting reads output](usage-recruit.md#output).

## Posterior TSV

Applies to: `process-mitochondrial-bam call`.

The tab-separated file contains a header and one row per read processed by the abundance model.

| Column | Description |
|:--------|:-------------|
| `read_name` | Read identifier. |
| `p_chrM` | Estimated probability that the read originated from mitochondrial DNA. |

Rows are sorted by read name to support reproducible comparisons.

Back to [Assigning mitochondrial reads output](usage-call.md#output).
