# Generate a panel of normals

`generate_pon` combines BAM feature tables from normal samples into a Panel of
Normals (PON) VCF.
Use the resulting indexed VCF with `--pon-vcf` when computing or filtering
variants with PON-derived features.

## How it works

For each normal sample, run `compute_bam_features panel-of-normals` to produce
a BAM features TSV.
List those TSV paths in a text file, one path per line, and pass the list to
`generate_pon`.

The command retains variants whose `duplex_af` is within the configured range
and whose `mapq_mean` meets the minimum threshold.
It aggregates retained variants across samples and writes their occurrence and
duplex allele-frequency statistics to the output VCF.

## Example command

Generate one feature table per normal sample:

```bash
compute_bam_features panel-of-normals \
    --genome REFERENCE.fa \
    --bam-input NORMAL_1.bam \
    --output-file NORMAL_1_BAM_FEATURES.tsv \
    --threads 16

compute_bam_features panel-of-normals \
    --genome REFERENCE.fa \
    --bam-input NORMAL_2.bam \
    --output-file NORMAL_2_BAM_FEATURES.tsv \
    --threads 16
```

Create the feature-list file and generate the PON:

```bash
printf '%s\n' \
    NORMAL_1_BAM_FEATURES.tsv \
    NORMAL_2_BAM_FEATURES.tsv > PON_FEATURES.txt

generate_pon \
    --feature-list PON_FEATURES.txt \
    --output-file pon.vcf.gz \
    --threads 16
```

## Input

| Input | Required | Notes |
|:---|:---|:---|
| Feature-list file | Yes | Passed with `--feature-list`. A non-empty text file containing one BAM features TSV path per line. |
| BAM features TSVs | Yes | Produced by `compute_bam_features panel-of-normals`. Each file must contain `chrom`, `pos`, `ref`, `alt`, `duplex_af`, and `mapq_mean` columns. |

## Output

| Output | Description |
|:---|:---|
| PON VCF | Written to `--output-file`. Contains one record per retained variant with occurrence and duplex allele-frequency INFO fields. See [Panel of Normals VCF](file-formats.md#panel-of-normals-vcf). |
| PON VCF index | Written next to a `.vcf.gz` output with the `.csi` suffix. |

Compressed output is BGZF-compressed and indexed automatically.

## Tips for analysis

Use normal samples that are representative of the sequencing protocol and
analysis workflow used for the samples that the PON will filter.
The duplex allele-frequency thresholds are inclusive.
Variants below `--min-duplex-af`, above `--max-duplex-af`, or below
`--min-mapq-mean` are excluded before aggregation.

## CLI options

Required parameters are shown in **bold**.

This command also accepts the common `-h, --help`, `-v, --version`, and
`-l, --log-level` options.

### Input options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--feature-list** | Path to a text file containing one BAM features TSV path per line. | Path to an existing non-empty text file |

### Output options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --output-file | Path to the output PON VCF. | Non-existent path [default: `pon.vcf.gz`] |

### PON filtering options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --min-duplex-af | Minimum `duplex_af`, inclusive, for retaining a variant. | Float from `0.0` to `1.0` [default: `0.01`] |
| --max-duplex-af | Maximum `duplex_af`, inclusive, for retaining a variant. Must be greater than or equal to `--min-duplex-af`. | Float from `0.0` to `1.0` [default: `0.1`] |
| --min-mapq-mean | Minimum `mapq_mean`, inclusive, for retaining a variant. | Float greater than or equal to `0.0` [default: `10.0`] |

### Performance options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --threads | Number of threads used (`0` = available hardware threads). The command uses no more threads than input feature files. | Integer greater than or equal to `0` [default: `1`] |

### Developer options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --warn-as-error | Treat warning messages as errors. | Boolean [default: `false`] |
