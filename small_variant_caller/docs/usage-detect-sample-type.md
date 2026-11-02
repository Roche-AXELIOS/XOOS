# Detect sample type

`detect_sample_type` infers whether a tumor-normal Mutect2 VCF corresponds to an
FFPE or a cell-line sample. It is the standalone entry point for the same
sample-type inference that `filter_variants tumor-normal-wgs` runs when
`--sample-type auto-detect` is used with automatic model selection.

## How it works

- The input VCF is read via [htslib](https://github.com/samtools/htslib).
- Low-allele-frequency tumor SNV observations (allele frequency at or below
  `--max-af`) are collected as error observations.
- A VCF phred-like metric is estimated from those observations over a
  denominator length. When `--target-bed` is provided, only records overlapping
  the BED intervals are evaluated and the BED interval length is used as the
  denominator. Otherwise, the denominator is the sum of `##contig` lengths in
  the VCF header, falling back to 3 Gbp when contig lengths are absent or sum to
  zero.
- Phred values at or below `--max-ffpe-phred` are classified as `ffpe`; values
  above the threshold are classified as `cell-line`.
- The inferred sample type is written to stdout, and the metrics behind the
  decision are logged.

## Example command

```bash
detect_sample_type \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --threads 16
```

The inferred sample type (`ffpe` or `cell-line`) is printed to stdout, which can
be captured to drive downstream model selection.

## Input

| Input | Required | Notes |
|:---|:---|:---|
| VCF file | Yes | Passed with `--vcf-input`. Tumor-normal Mutect2 VCF file. Index required for parallel scanning. |
| Target regions | No | Passed with `--target-bed`. Limits detection to the BED intervals. |

## Output

| Output | Description |
|:---|:---|
| Sample type | Written to stdout: `ffpe` or `cell-line`. |

## CLI options

Required parameters are shown in **bold**.

This command also accepts the common `-h, --help`, `-v, --version`, and
`-l, --log-level` options. Unlike the other SVC submodules, it does not accept a
`--config` file.

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--vcf-input** | Path to input tumor-normal VCF file. | Path to an existing VCF file |
| --target-bed | Limits sample-type detection to 0-based target regions. | An existing BED file (uncompressed or gzip/BGZF-compressed `.gz`) |
| --max-af | Maximum allele frequency for SNV error observations used by the detector. | Fraction [default: `0.10`] |
| --max-ffpe-phred | Maximum VCF phred classified as FFPE. Values above this threshold are classified as cell-line. | Non-negative float [default: `42.5`] |
| --threads | Number of threads used (`0` = available hardware threads). Indexed VCF required for parallel scanning. | Integer >= 0 [default: `1`] |
| --warn-as-error | Treat warning messages as errors. | [default: `false`] |
