# Assigning mitochondrial reads

## How it works

`process-mitochondrial-bam call` labels GAF placements as mitochondrial or nuclear by comparing their graph nodes with a mitochondrial path-node set.
It initializes a 2-class abundance model from mitochondrial and autosomal coverage, applies length normalization, and estimates a mitochondrial posterior probability for each read.

The command retains mitochondrial placements for assigned reads and replaces GAF column 12 with a capped Phred quality score derived from the posterior probability.
The end-to-end workflow surjects this GAF file after assignment so origin-spanning mitochondrial alignments receive consistent supplementary records and `SA:Z` tags.

## Example command

```shell
process-mitochondrial-bam call \
  --gaf all.gaf \
  --paths-gaf resources/chrm_paths.grch38.gaf \
  --orig-bam sample.bam \
  --recruit-regions recruited.bed \
  --nr-positions resources/numt_catalog.grch38.nr.tsv \
  --control-regions resources/control_regions.grch38.bed \
  --out-gaf chrM.gaf \
  --threads 8
```

## Input

| Input | Required | Notes |
|:-------|:---------|:------|
| Combined GAF | Yes | Multimapping graph alignments for mitochondrial and recruited reads. |
| Paths GAF | Yes | GAF line that defines the mitochondrial graph-path node set. |
| Original BAM | Yes | Coordinate-sorted, indexed, NR-tagged BAM file used for coverage initialization and polymorphic NuMT presence checks. |
| Recruited regions | No | Optional merged BED file used for non-mitochondrial length normalization. |
| NuMT anchors | No | Optional NuMT anchor TSV used for presence-checked insertion lengths. |
| Control regions | No | Optional autosomal BED file used to estimate background depth. |

## Output

The primary output is a GAF file containing retained mitochondrial placements with updated mapping qualities.
Use `vg surject`, followed by BAM sorting and indexing, to create a final mitochondrial BAM file.

The optional [posterior TSV](file-formats.md#posterior-tsv) reports each read's estimated mitochondrial probability.

## Tips for analysis

Ensure that `--paths-gaf` is consistent with the GBZ used to produce `--gaf`.
Supply `--recruit-regions` to activate length normalization.
Depth overrides bypass the corresponding coverage estimate and support controlled evaluation.

## CLI options

This command also accepts the common `-h, --help`, `-v, --version`, and `-l, --log-level` options.

### Input options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --gaf | Combined multimapping GAF file. | File path; required |
| --paths-gaf | Paths GAF file that defines the mitochondrial node set. | File path; required |
| --orig-bam | Original coordinate-sorted, indexed BAM file. | File path; required |
| --recruit-regions | Recruited-region BED file for length normalization. | File path; optional |
| --nr-positions | [NuMT anchor TSV](file-formats.md#numt-anchor-tsv) with insertion lengths. | File path; optional |
| --control-regions | Autosomal control-region BED file. | File path; optional |

### Output options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --out-gaf | Mitochondrial GAF file with posterior-derived mapping qualities. | File path; required |
| --posteriors-out | [Posterior TSV](file-formats.md#posterior-tsv) for evaluation and debugging. | File path; optional |

### Assignment options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --read-len | Read length used for capture-band padding. | Integer; optional; default: estimated from GAF |
| --window | Presence-check radius in bp around each catalog anchor. | Integer; optional; default: `10` |
| --min-mapq | Starting per-base mapping-quality threshold for mitochondrial depth. | Integer; optional; default: `55` |
| --chrm-depth | Mitochondrial depth override. | Value; optional; default: estimated from BAM |
| --non-chrm-depth | Non-mitochondrial depth override. | Value; optional; default: estimated from control regions |
| --chrm-len | Mitochondrial contig length. | Integer; optional; default: BAM header |
| --max-iter | Maximum number of expectation-maximization iterations. | Integer; optional; default: `100` |
| --tol | Expectation-maximization convergence tolerance. | Value; optional; default: `1e-06` |
| --chrm-name | Mitochondrial graph-path and contig name. | Name; optional; default: `chrM` |

### Performance options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --threads | Number of threads; `0` uses the available hardware threads. | Integer; optional; default: `1` |
