# Usage

Use this page to choose the right Alignment Metrics command.
Follow the links below to navigate to the full command guides.

## Choosing between `all` and the specific commands

Use `all` when you want the normal Alignment Metrics deliverable and one shared set of inputs, filters, and output directories.
Use a specific command when you need only one metrics family,
when a reference is unavailable for a pure `read` run,
or when you want a faster debugging loop on one output family.

## Compute all metrics

Run the complete Alignment Metrics workflow.
Input: **BAM + reference + optional BED** -> Output: **`coverage/`, `accuracy/`, and `read_metrics/` directories**.

```shell
alignment_metrics all --bam INPUT.bam --reference REFERENCE.fa --output-dir OUTPUT_DIR
```

To run the complete workflow, go to [Compute all metrics](usage-all.md).

For label-stratified runs with `all`,
build the labeled BED input with [`merge_and_label_bed`](usage-merge-and-label-bed.md).

## Compute coverage metrics

Compute coverage and optional targeted-enrichment or homopolymer coverage metrics.
Input: **BAM + reference + optional BED** -> Output: **`coverage/` directory**.

```shell
alignment_metrics coverage --bam INPUT.bam --reference REFERENCE.fa --output-dir OUTPUT_DIR
```

To run coverage analysis, go to [Compute coverage metrics](usage-coverage.md).

For label-stratified runs with `coverage`,
build the labeled BED input with [`merge_and_label_bed`](usage-merge-and-label-bed.md).

## Compute accuracy metrics

Compute base-level accuracy and error-breakdown metrics.
Input: **BAM + reference + optional BED** -> Output: **`accuracy/` directory**.

```shell
alignment_metrics accuracy --bam INPUT.bam --reference REFERENCE.fa --output-dir OUTPUT_DIR
```

To run accuracy analysis, go to [Compute accuracy metrics](usage-accuracy.md).

For label-stratified runs with `accuracy`,
build the labeled BED input with [`merge_and_label_bed`](usage-merge-and-label-bed.md).

## Compute GC bias metrics

Compute GC bias metrics.
Input: **BAM + reference + optional BED** -> Output: **`gc_bias/` directory**.

```shell
alignment_metrics gcbias --bam INPUT.bam --reference REFERENCE.fa --output-dir OUTPUT_DIR
```

To run GC bias analysis, go to [Compute GC bias metrics](usage-gcbias.md).

For label-stratified runs with `gcbias`,
build the labeled BED input with [`merge_and_label_bed`](usage-merge-and-label-bed.md).

## Compute read metrics

Compute read-count and read-length metrics without requiring a reference.
Input: **BAM + optional BED** -> Output: **`read_metrics/` directory**.

```shell
alignment_metrics read --bam INPUT.bam --output-dir OUTPUT_DIR
```

To run read metrics analysis, go to [Compute read metrics](usage-read.md).

For label-stratified runs with `read`,
build the labeled BED input with [`merge_and_label_bed`](usage-merge-and-label-bed.md).
