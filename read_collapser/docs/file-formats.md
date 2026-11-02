# File formats

This page is the reference for the tool-specific formats that Read Collapser
produces. Usage details for `markdup` and `consensus` are on the
[Mark duplicates](usage-markdup.md) and [Consensus](usage-consensus.md) pages,
and metric file schemas are on the [Metrics](metrics.md) page.

## Custom BAM tags

Applies to: `markdup`, `consensus` with `--write-cluster-bam`.

When `markdup` assigns a read to a cluster (or when `consensus` is run with `--write-cluster-bam`), each read receives two auxiliary tags:

| Tag     | Type        | Definition                                                                                                     |
|:--------|:------------|:---------------------------------------------------------------------------------------------------------------|
| DI      | string `Z`  | A unique cluster ID. The numerical value is arbitrary and only serves to distinguish one cluster from another. |
| DS      | integer `i` | The number of reads in the cluster.                                                                            |

In `markdup`, these tags can be disabled with `--exclude-cluster-tags`.
In `consensus`, they are only written when `--write-cluster-bam` is enabled; duplicate flags (0x400) from the input BAM are cleared in the cluster BAM output, and the cluster BAM is not sorted or indexed.

Back to [Mark duplicates output](usage-markdup.md#output) and
[Consensus output](usage-consensus.md#output).

## Consensus FASTQ record format

Applies to: `consensus`.

Each consensus read is represented in the FASTQ output file as follows:

```text
@read_name
consensus_sequence (for example, CACTGA)
+
base_qualities (for example, III!II)
```

The `read_name` follows the format
`[cluster_id]-[num_partial_fwd]-[num_partial_rev]-[num_full_fwd]-[num_full_rev]-[effective_cluster_size]`,
where `-` is the delimiter. The `base_qualities` are ASCII encoded Phred quality
scores.

Cluster size vs effective cluster size:

| Concept | Definition | Used for |
|:--------|:-----------|:---------|
| Cluster size | Count of reads assigned to a cluster during clustering (pre-consensus). | Cluster size distributions and filtering (`--min-cluster-size`, `--max-cluster-size`). |
| Effective cluster size | Average read support per consensus position (total bases in consensus matrix / covered positions). This can be lower than the raw cluster size when coverage is uneven. | Consensus read name annotation and strand-aware filtering (`--min-same-strand-cluster-size`, `--min-mixed-strand-cluster-size`). |

Why this distinction matters: cluster size tells you how many reads were
grouped, while effective cluster size tells you how much support each consensus
base actually had. Using both lets you tune clustering and consensus quality
independently and interpret metrics correctly.

| Read name component | Definition |
|:--------------------|:-----------|
| cluster_id | The consensus read's unique cluster ID. The numerical value is arbitrary and only serves to distinguish one cluster from another. |
| num_partial_fwd | The number of partial forward reads in the cluster used to generate the consensus read. |
| num_partial_rev | The number of partial reverse reads in the cluster used to generate the consensus read. |
| num_full_fwd | The number of full forward reads in the cluster used to generate the consensus read. |
| num_full_rev | The number of full reverse reads in the cluster used to generate the consensus read. |
| effective_cluster_size | Average read support per consensus position (total bases in consensus matrix / covered positions). |

Because Read Collapser cannot determine direction for unmapped reads, the
`read_name` always follows the format
`cluster_id-0-0-0-0-effective_cluster_size`, where
`effective_cluster_size` is 1 unless `--duplex-library-type` is set to
`parent-parent` or `parent-daughter`, in which case `effective_cluster_size` is
approximated with
`1 + round(concordant bases / (concordant bases + simplex bases + discordant bases))`.

All consensus reads are reported as forward strand (`+`) because Read Collapser
assumes the aligner that mapped the reads in the input BAM took the
reverse complement of any reads mapping to the reverse strand. Under this
assumption, Read Collapser only considers sequences based on the forward strand.

Back to [Consensus output](usage-consensus.md#output).

## Consensus auxiliary tags

Applies to: `consensus` when `--include-consensus-aux-tags` is enabled.

These duplex consensus auxiliary tags let downstream pipelines built for this tag schema consume Read Collapser output.
They provide per-position and summary insight into how many raw reads contributed to each consensus base and how many
disagreed with the final call.

`a` refers to the forward single-strand consensus, `b` to the reverse single-strand consensus, and `c` to the combined
(both-strand) consensus. The tags are appended to the read-name line of each consensus FASTQ record, separated from the
read name and from each other by tab characters (`\t`), using standard SAM auxiliary field encoding.

Per-position arrays (length equal to the consensus sequence length):

| Tag | Type      | Definition                                                                                       |
|:-----|:-----------|:--------------------------------------------------------------------------------------------------|
| ad  | array `B:s` | The count of forward raw-read bases contributing to the consensus at each position.            |
| bd  | array `B:s` | The count of reverse raw-read bases contributing to the consensus at each position.            |
| ae  | array `B:s` | The count of forward raw-read bases disagreeing with the consensus base at each position.      |
| be  | array `B:s` | The count of reverse raw-read bases disagreeing with the consensus base at each position.      |

The per-position arrays are only emitted for strands present in the cluster. A forward-only cluster omits the reverse
arrays (`bd`/`be`), and a reverse-only cluster omits the forward arrays (`ad`/`ae`), because the absent strand's
per-position counts would be uniformly zero. Mixed-strand clusters emit all four arrays. The per-strand and combined
summary tags below are always emitted regardless of cluster strandedness.

Per-strand and combined summaries:

| Tag        | Type       | Definition                                                                                      |
|:------------|:------------|:-------------------------------------------------------------------------------------------------|
| aD, bD, cD | integer `i` | The maximum raw-read depth at any position in the (forward / reverse / combined) consensus.     |
| aM, bM, cM | integer `i` | The minimum raw-read depth at any covered position in the (forward / reverse / combined) consensus. |
| aE, bE, cE | float `f`   | The fraction of raw-read bases disagreeing with the final (forward / reverse / combined) consensus calls. |

Notes:

- Due to the type of the per-base arrays being `B:s` (array of signed 16-bit integers), the maximum value for any position is `32767`.
  Any counts above this value are capped at `32767`.
- The `ac`, `bc` (single-strand consensus bases) and `aq`, `bq` (single-strand consensus qualities) tags are not
  emitted, because Read Collapser computes a single combined consensus rather than two independent single-strand
  consensuses, so per-strand base and quality strings are not available.
- The minimum-depth tags (`aM`/`bM`/`cM`) ignore positions with zero depth on that strand, so a strand that does not
  cover part of the consensus does not force the minimum to zero.
- The error rate is `total disagreeing bases / total depth` over all positions on that strand (or combined). It is `0`
  when the strand has no coverage.

The read name with these tags follows this format (shown for a mixed-strand cluster, which emits all four per-position
arrays):

```text
[cluster_id]-[num_partial_fwd]-[num_partial_rev]-[num_full_fwd]-[num_full_rev]-[effective_cluster_size]   ad:B:s,...   bd:B:s,...   ae:B:s,...   be:B:s,...   aD:i:...   bD:i:...   cD:i:...   aM:i:...   bM:i:...   cM:i:...   aE:f:...   bE:f:...   cE:f:...
```

A forward-only cluster omits `bd`/`be`, and a reverse-only cluster omits `ad`/`ae`; the summary tags
(`aD`/`bD`/`cD`, `aM`/`bM`/`cM`, `aE`/`bE`/`cE`) are always present.

Back to [Consensus output](usage-consensus.md#output).

The metric files (`metrics/summary_stats.tsv`, `metrics/cluster_size_distribution_summary.tsv`,
and `metrics/cluster_size_distributions.tsv`) are documented on the [Metrics](metrics.md) page.
