# Running the end-to-end workflow

## How it works

The `run-process-mitochondrial-bam.sh` helper script reproduces the chrM realignment portion of the `process-mitochondrial-bam` workflow as implemented in the [xoosnf pipeline](../../docs/xoosnf/README.md) and documented [here](../../docs/xoosnf/mitochondrial-bam-post-processing.md). This helper script can be found within the `bam_post_process` image. **NOTE**: The vg version in this image may not match that of the xoosnf pipeline.
It performs the following stages:

1. `process-mitochondrial-bam recruit` selects mitochondrial and NuMT-associated alignments and writes a FASTQ file.
2. `vg giraffe` realigns the selected reads, retaining multiple candidate placements.
3. `process-mitochondrial-bam call` estimates mitochondrial origin probabilities and rewrites GAF mapping qualities.
4. `vg surject` converts retained graph alignments to BAM records, and `samtools` sorts and indexes the output.

Temporary files remain in a temporary working directory and are removed when the script exits.

## Example command

```shell
run-process-mitochondrial-bam.sh \
  --orig-bam sample.bam \
  --gbz graph.gbz \
  --min graph.min \
  --dist graph.dist \
  --zipcodes graph.zipcodes \
  --out-bam sample.chrM.bam \
  --threads 8
```

## Input

| Input | Required | Notes |
|:-------|:---------|:------|
| Original BAM | Yes | Coordinate-sorted, indexed, NR-tagged BAM file produced by pangenome alignment and surjection. |
| GBZ | Yes | Pangenome graph that matches the original alignment and companion indexes. |
| Minimizer index | Yes | Giraffe `.min` index for the GBZ. |
| Distance index | Yes | Giraffe `.dist` index for the GBZ. |
| Zipcode index | No | Optional `.zipcodes` index that matches the minimizer index. |

The image supplies GRCh38 NuMT, control-region, and mitochondrial-path resources by default.

## Output

The workflow writes a coordinate-sorted mitochondrial BAM file to `--out-bam` and an adjacent BAI index.
The output includes reads assigned to `chrM`, with mapping qualities derived from their posterior mitochondrial probabilities.

## Tips for analysis

Use graph indexes from the same graph build.
The packaged mitochondrial-path resource is specific to its documented GRCh38 pangenome graph.

## CLI options

The script accepts `-h, --help` in addition to the following options.

### Input options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --orig-bam | Original coordinate-sorted, indexed, NR-tagged BAM file. | File path; required |
| --gbz | Pangenome graph in GBZ format. | File path; required |
| --min | Giraffe minimizer index. | File path; required |
| --dist | Giraffe distance index. | File path; required |
| --zipcodes | Giraffe zipcode index. | File path; optional |
| --nr-positions | NuMT anchor TSV file. | File path; optional; default: packaged resource |
| --ucsc-bed | Reference NuMT BED file. | File path; optional; default: packaged resource |
| --control-bed | Autosomal control-region BED file. | File path; optional; default: packaged resource |

### Output options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --out-bam | Output mitochondrial BAM file. | File path; required |

### Assignment options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --anchor-pad | Padding in bp around each polymorphic NuMT anchor. | Integer; optional; default: `150` |
| --chrm-name | Mitochondrial contig and graph-path name. | Name; optional; default: `chrM` |
| --window | Proximity window in bp around NuMT anchors. | Integer; optional; default: `150` |
| --multimap | Maximum number of placements requested from `vg giraffe`. | Integer; optional; default: `10` |

### Performance options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --threads | Number of threads used by the workflow. | Integer; optional; default: `4` |
