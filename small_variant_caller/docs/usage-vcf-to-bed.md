# Convert VCF to BED

`vcf_to_bed` builds a regions BED file from the variant positions in a VCF file.
It is a utility helper: use it to derive the target regions that speed up BAM
feature computation for model training with
[compute_bam_features](usage-compute-bam-features.md).

## How it works

- The input VCF is read via [htslib](https://github.com/samtools/htslib).
- Each variant contributes an interval derived from its position, optionally
  extended by `--left-pad` and `--right-pad`.
- Intervals within `--collapse-distance` of each other are merged after padding
  is applied.
- The result is written as a BED file of 0-based regions.

## Example command

```bash
vcf_to_bed \
    --vcf-input INPUT.vcf.gz \
    --output-file VCF_REGIONS.bed \
    --left-pad 0 \
    --right-pad 0 \
    --collapse-distance 0
```

The output BED can be passed as `--target-regions` to
[compute_bam_features](usage-compute-bam-features.md), restricting feature
computation to the variant loci.

## Input

| Input | Required | Notes |
|:---|:---|:---|
| VCF file | Yes | Passed with `--vcf-input`. The input VCF file. |
| Target regions | No | Passed with `--target-regions`. Optional BED to restrict output regions. |

## Output

| Output | Description |
|:---|:---|
| Regions BED | Written to `--output-file`. BED file of 0-based variant regions. |

## CLI options

Required parameters are shown in **bold**.

This command also accepts the common `-h, --help`, `-v, --version`, and
`-l, --log-level` options.

### Input options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--vcf-input** | Path to the input VCF file. | Path to an existing VCF file |

### Output options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --output-file | Output BED file name. | [default: `vcf-regions.bed`] |

#### Output BED shaping options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --left-pad | Left-padding to variant start position for the output BED file. | Integer >= 0 [default: `0`] |
| --right-pad | Right-padding to variant end position for the output BED file. | Integer >= 0 [default: `0`] |
| --collapse-distance | Distance to collapse nearby variants after padding is applied for the output BED file. | Integer >= 0 [default: `0`] |

### Region options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --target-regions | Restricts output to 0-based target regions. | An existing BED file (uncompressed or gzip/BGZF-compressed `.gz`) |

### Performance options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --threads | Number of threads used (`0` = available hardware threads). | Integer >= 0 [default: `1`] |

### Developer options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --warn-as-error | Treat warning messages as errors. | [default: `false`] |
