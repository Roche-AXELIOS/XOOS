# Usage

Read Collapser provides two main subcommands, `markdup` and `consensus`, plus a
`read_collapser_rescue_secondary` command for pangenome-aligned data. This page
summarizes each command and links to its dedicated Usage page, where you will
find how it works, an example command, input and output details, tips, and the
full CLI options.

## Mark duplicates

Identifies, marks, and optionally removes duplicate reads. Within each cluster,
all reads are treated as duplicates except the single highest-quality
representative. Input: **position-sorted BAM** → Output: **duplicate-marked
BAM**.

```bash
read_collapser markdup --preset wgs-duplex --bam HG001.bam
```

For inputs, outputs, tips, and the full CLI options, go to the
[Mark duplicates](usage-markdup.md) page.

{% hint style="warning" %}
For pangenome-aligned data, run
[`read_collapser_rescue_secondary`](usage-rescue-secondary.md) before
duplicate marking to recover reads whose primary alignment could not be projected
to the linear reference.
{% endhint %}

For inputs, outputs, tips, and the full CLI options, go to the
[Rescue secondary alignments](usage-rescue-secondary.md) page.

## Consensus

Generates high-quality consensus sequences by combining evidence from all reads
in a cluster through a matrix-based majority vote, correcting random sequencing
errors and damage artifacts. Input: **position-sorted BAM** → Output:
**consensus FASTQ**.

```bash
read_collapser consensus --preset wgs-duplex-cfdna --bam HG001.bam
```

For inputs, outputs, tips, and the full CLI options, go to the
[Consensus](usage-consensus.md) page.
