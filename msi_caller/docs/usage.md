# Usage

Use this page to choose an MSI calling workflow or generate a panel-of-normals background model.

The `msi_caller` application provides the `tumor-only` and `tumor-normal` subcommands.
The separate `background_generator` application creates the panel-of-normals model required by `tumor-only` analysis.

The `--threads` and `--min-msi-classification-fraction` options are optional, but recommended for caller commands.
The examples below use `--threads 8` and the recommended `--min-msi-classification-fraction 0.3`.

## Tumor-only analysis

Compare tumor repeat-length distributions with a panel-of-normals background model.
Input: **indexed tumor BAM + background model JSON** → Output: **VCF calls + JSON report + TSV metrics**.

```shell
msi_caller tumor-only \
  --tumor-bam TUMOR.bam \
  --background-threshold-json BACKGROUND.json \
  --output-dir OUTPUT_DIR \
  --min-msi-classification-fraction 0.3 \
  --threads 8
```

To run tumor-only analysis, go to [Tumor-only analysis](usage-tumor-only.md).

## Tumor-normal analysis

Compare repeat-length distributions from matched tumor and normal BAMs.
Input: **indexed tumor and normal BAMs + sites TSV** → Output: **VCF calls + JSON report + TSV metrics**.

```shell
msi_caller tumor-normal \
  --tumor-bam TUMOR.bam \
  --normal-bam NORMAL.bam \
  --microsatellite-sites-thresholds-tsv SITES.tsv \
  --output-dir OUTPUT_DIR \
  --min-msi-classification-fraction 0.3 \
  --threads 8
```

To run tumor-normal analysis, go to [Tumor-normal analysis](usage-tumor-normal.md).

## Generate a background model

Create a panel-of-normals background model for tumor-only analysis.
Input: **indexed BAM list + sites TSV** → Output: **background model JSON**.

```shell
background_generator \
  --pon-bam-list PON_BAMS.txt \
  --microsatellite-sites-thresholds-tsv SITES.tsv \
  --output-background-threshold-json BACKGROUND.json
```

To generate a background model, go to [Generate a background model](usage-background-generator.md).
