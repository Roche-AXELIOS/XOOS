# Convert features to VCF

`features_to_vcf` converts a BAM features TSV into a coordinate-sorted, indexed,
single-sample VCF. It is a utility helper: use it to turn the variant
coordinates produced by [compute_bam_features](usage-compute-bam-features.md)
into a FORMAT-bearing VCF that downstream tools requiring a VCF can read.

## How it works

- The input TSV must contain `chrom`, `pos`, `ref`, and `alt` columns; all other
  columns are ignored.
- Variants are sorted by reference FASTA contig order, then by position, then by
  alternate allele.
- Variants on contigs absent from the reference FASTA index are skipped with a
  warning.
- Each record is written with `FILTER=PASS` and synthetic diploid heterozygous
  genotype fields (`GT=0/1`, `AC=1`, `AF=0.5`, `AN=2`, `DP=2`, `AD=1,1`, and
  `GQ=50`) so that tools requiring a FORMAT-bearing VCF can read the variant
  sites. These fields are placeholders and do not represent measured genotype
  evidence.

## Example command

```bash
features_to_vcf \
    --input-file BAM_FEATURES.tsv \
    --genome REFERENCE.fa \
    --sample-name SAMPLE \
    --output-file features.vcf.gz
```

## Input

| Input | Required | Notes |
|:---|:---|:---|
| BAM features TSV | Yes | Passed with `--input-file`. Must contain `chrom`, `pos`, `ref`, and `alt` columns. |
| Reference genome | Yes | Passed with `--genome`. Indexed FASTA; provides contig order and filters unknown contigs. |

## Output

| Output | Description |
|:---|:---|
| Single-sample VCF | Written to `--output-file`. Coordinate-sorted, indexed VCF with synthetic heterozygous genotype fields. |

## CLI options

Required parameters are shown in **bold**.

This command also accepts the common `-h, --help`, `-v, --version`, and
`-l, --log-level` options.

### Input options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--input-file** | Path to the input BAM features TSV file. | Path to an existing non-empty file |
| **--genome** | Path to the indexed reference genome FASTA. | Path to an existing FASTA file |

### Output options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --output-file | Output VCF file name. | [default: `features.vcf.gz`] |
| --sample-name | Sample name to write in the output VCF header. | String [default: `SAMPLE`] |

### Developer options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --warn-as-error | Treat warning messages as errors. | [default: `false`] |
