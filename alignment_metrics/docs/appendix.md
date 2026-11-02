# Appendix

## Concepts and terminology

### General definitions

| Term | Definition |
|:------|:------------|
| Read | A single sequencing read from the BAM file. |
| Alignment | A single alignment of a read to the reference genome. A read can have multiple alignments, including primary, secondary, and supplementary alignments. |
| Mapped read | A read that has at least one alignment to the reference genome. |
| Unmapped read | A read that is not aligned to the reference genome. This includes placed and unplaced unmapped reads. |
| Forward read | A read aligned in the forward direction relative to the reference genome. |
| Reverse read | A read aligned in the reverse direction relative to the reference genome. |
| Full read | A read with a complete insert. Alignment Metrics currently uses the presence of both expected UMIs as a proxy. |
| Partial read | A read with an incomplete insert. Alignment Metrics currently uses the absence of one expected UMI as a proxy. |

### Definitions for consensus reads

| Term | Definition |
|:------|:------------|
| Cluster size | Number of reads in the cluster from which the consensus read is generated. |
| Forward cluster read | Consensus read produced from a cluster containing only forward reads. |
| Reverse cluster read | Consensus read produced from a cluster containing only reverse reads. |
| Mixed-strand cluster read | Consensus read produced from a cluster containing both forward and reverse reads. |
| Full cluster read | Consensus read produced from a cluster containing only full reads. |
| Mixed-full-and-partial cluster read | Consensus read produced from a cluster containing both full and partial reads. |
| Partial cluster read | Consensus read produced from a cluster containing only partial reads. |

### Definitions for duplex reads

| Term | Definition |
|:------|:------------|
| R1 | The region of a duplex read between the start adapter and the mid adapter. |
| R2 | The region of a duplex read between the mid adapter and the end adapter. |
| Concordant duplex base | Base that agrees between R1 and R2. |
| Discordant duplex base | Base that disagrees between R1 and R2. |
| Simplex base | Base present in only one of R1 or R2. |
| YC tag | BAM tag that encodes the two constituent reads of a duplex consensus read. |

Structure of a duplex read:

![SBX-D adapter structure](assets/sbx-d-adapter-text-outlined.svg)

Example of simplex and discordant positions:

```text
  simplex bases
    ┌─┴──┐
R1: AAGTCGACGTACTTGAAATGT
          |||||x|||||||||
R2:       ACGTATTTGAAATGT
               ↑
        discordant base
```

### Definitions for homopolymer accuracy

| Term | Definition |
|:------|:------------|
| Homopolymer | Stretch of repeating identical bases of length greater than or equal to 2. |
| Homopolymer region | Region in the reference genome that contains a homopolymer. |
| Spanning read | Read that covers the homopolymer region with at least one extra flanking base on both sides. |
| Effective read | Spanning read where all bases aligned to the homopolymer region pass the quality and base-type filters. Unless otherwise specified, insertions must also be homogeneous. |
| Homogeneous insertion | Insertion in a homopolymer region of the same base as the homopolymer base. |
| Heterogeneous insertion | Insertion in a homopolymer region of a different base from the homopolymer base. |

## Performance and cost notes

### `all`

- Runtime scales with the size of the queried regions, depth of coverage, and whether homopolymer metrics are enabled.
- `--enable-target-enrichment-metrics` disables region partitioning to preserve target-level summaries, which can increase peak memory for very large BED intervals.
- Decreasing `--region-size` reduces per-thread memory demand, but too small a region size can increase overhead.

### `coverage`

- Coverage-only runs are typically faster than `all` because they do not compute read or accuracy outputs.
- Enabling homopolymer metrics adds reference traversal and homopolymer bookkeeping.
- Enabling target-enrichment metrics disables region partitioning, which can increase peak memory for very large intervals.

### `accuracy`

- Accuracy runs require the reference and generally cost more than read-only runs.
- Homopolymer metrics increase runtime and memory because the run tracks homopolymer contexts and effective-read classifications.
- Smaller `--region-size` values reduce per-thread memory demand at the cost of extra scheduling overhead.

### `read`

- `read` is the cheapest Alignment Metrics command because it does not need the reference and does not build the full accuracy stack.
- `--enable-target-enrichment-metrics` can increase peak memory for very large BED intervals because region partitioning is disabled for target-level summaries.
- Read trimming affects both the counts and the read-length histograms.

### `merge_and_label_bed`

- The helper is typically cheap compared with the main Alignment Metrics run.
- The primary practical limit is the label surface. At most 16 distinct labels are supported.
