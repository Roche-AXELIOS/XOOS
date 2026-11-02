# Recruiting reads

## How it works

`process-mitochondrial-bam recruit` selects reads through 2 static NuMT catalog routes.
The interval route selects reads that overlap padded polymorphic NuMT anchors or reference NuMT intervals.
The NR-proximity route selects reads whose `NR` tag indicates the position of a polymorphic NuMT anchor.

The command can write the selected names, merged catalog intervals, and single-end FASTQ records during one scan of the original BAM file.

## Example command

```shell
process-mitochondrial-bam recruit \
  --orig-bam sample.bam \
  --nr-positions resources/numt_catalog.grch38.nr.tsv \
  --ucsc-bed resources/nuMtSeq_hg38.ucsc.bed \
  --out-names recruited.names \
  --out-regions recruited.bed \
  --out-fastq recruited.fastq \
  --threads 8
```

## Input

| Input | Required | Notes |
|:-------|:---------|:------|
| Original BAM | Yes | Coordinate-sorted BAM file with an adjacent BAI index and `NR:Z` tags on off-reference alignments. |
| NuMT anchors | No | Tool-specific TSV containing polymorphic NuMT anchors. |
| Reference NuMT regions | No | BED file containing reference NuMT intervals. |

Use `--catalog-regions` to replace the derived interval set with a pre-merged BED file.

## Output

The command can produce the following outputs:

- A sorted, deduplicated [read-name list](file-formats.md#recruited-read-name-list).
- A merged BED file for `process-mitochondrial-bam call --recruit-regions`.
- A single-end FASTQ file containing mitochondrial and recruited primary reads.

If `--out-names` is omitted, the read-name list is written to standard output.

## CLI options

This command also accepts the common `-h, --help`, `-v, --version`, and `-l, --log-level` options.

### Input options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --orig-bam | Original BAM file. | File path; required |
| --catalog-regions | Pre-merged NuMT BED file that replaces interval derivation. | File path; optional |
| --nr-positions | [NuMT anchor TSV](file-formats.md#numt-anchor-tsv) for interval and NR-proximity recruitment. | File path; optional |
| --ucsc-bed | Reference NuMT BED file for interval recruitment. | File path; optional |

### Output options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --out-names | Recruited read-name list. | File path; optional; default: standard output |
| --out-regions | Merged NuMT interval BED file. | File path; optional |
| --out-fastq | Single-end FASTQ file for graph realignment. | File path; optional |

### Recruitment options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --anchor-pad | Padding in bp around each NuMT anchor for interval recruitment. | Integer; optional; default: `150` |
| --chrm-name | Mitochondrial contig name used when writing FASTQ records. | Name; optional; default: `chrM` |
| --window | Proximity window in bp around catalog NR positions. | Integer; optional; default: `150` |

### Performance options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --threads | Number of threads; `0` uses the available hardware threads. | Integer; optional; default: `1` |
