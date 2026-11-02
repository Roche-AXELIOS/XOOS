# Usage

The Small Variant Caller provides nine commands: `filter_variants`,
`phase_variants`, `train_model`, `compute_bam_features`, `compute_vcf_features`,
`generate_pon`, `vcf_to_bed`, `detect_sample_type`, and `features_to_vcf`.
This page summarizes each command and links to its dedicated Usage page, where
you will find how it works, an example command, input and output details, tips,
and the full CLI options.

For an overview of the tool, example applications, and system requirements, see
the [Small Variant Caller overview](../README.md#overview).

**Before you begin (for `filter_variants`, `phase_variants`, `compute_bam_features`, and `compute_vcf_features`):**
Generate the required candidate VCF and realigned BAM by following
[Generating candidate calls](generating-candidate-calls.md).

## Filter variants

Applies pre-trained ML models to filter and re-genotype an upstream call set.
Input: **candidate VCF + BAM + reference** → Output: **filtered VCF**.

```bash
filter_variants germline-wgs \
    --bam-input GATK_OUTPUT.bam \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --genome REFERENCE.fa \
    --aligner bwa \
    --pop-af-vcf POPAF.vcf.gz \
    --output-dir /PATH/TO/OUTPUT/DIRECTORY \
    --vcf-output FILTER_VARIANTS.vcf.gz \
    --threads 16
```

For inputs, outputs, tips, and the full CLI options, see the
[Filter variants](usage-filter-variants.md).

For `filter_variants tumor-normal-wgs`, specify the tumor read-group with
`--tumor-sample-name TUMOR`. The option identifies tumor reads in the matched
tumor-normal BAM; no separate normal sample-name option is required.

## Phase variants

Uses read-level evidence from an aligned BAM to phase nearby variants and
optionally merge phased variants into MNV or complex records.
Input: **candidate VCF + indexed BAM** → Output: **phased VCF**.

```bash
phase_variants \
    tumor-only-wgs \
    --bam-input MUTECT2_BAMOUT.bam \
    --vcf-input SVC_OUTPUT.vcf.gz \
    --output-dir phased \
    --vcf-output PHASED_OUTPUT.vcf.gz \
    --tumor-sample-name TUMOR
```

For inputs, outputs, tips, and the full CLI options, see the
[Phase variants](usage-phase-variants.md).

## Train models

Trains LightGBM models for variant filtering and re-genotyping from
pre-computed feature files and a truth VCF.
Input: **BAM or VCF feature files + truth VCF** → Output: **model files**.
The supported germline WGS library types are SBX-Simplex (SBX-S), SBX-Duplex
(SBX-D), and SBX-Fast.
For the germline WGS `--library-type` option, see
[Library types for germline WGS](appendix.md#library-types-for-germline-wgs).

```bash
train_model germline-wgs \
    --library-type sbx-d \
    --positive-bam-features BAM_FEATURES_FILE_LIST.txt \
    --positive-vcf-features VCF_FEATURES_FILE_LIST.txt \
    --truth-vcfs TRUTH_VCF_FILE_LIST.txt \
    --snv-output-file SNV_MODEL.txt.gz \
    --indel-output-file INDEL_MODEL.txt.gz \
    --threads 16
```

For inputs, outputs, tips, and the full CLI options, see the
[Train models](usage-train-model.md).

## Compute BAM features

Computes per-variant ML features from one or more aligned BAM files.
Input: **BAM + reference (+ regions BED)** → Output: **BAM features TSV**.

```bash
compute_bam_features germline-wgs \
    --genome REFERENCE.fa \
    --bam-input GATK_OUTPUT.bam \
    --target-regions VCF_REGIONS.bed \
    --output-file BAM_FEATURES.tsv \
    --threads 16
```

For inputs, outputs, and the full CLI options, see the
[Compute BAM features](usage-compute-bam-features.md).

## Compute VCF features

Computes per-variant ML features from a candidate VCF.
Input: **VCF + reference (+ regions BED)** → Output: **VCF features TSV (+ regions BED)**.
The supported germline WGS library types are SBX-Simplex (SBX-S), SBX-Duplex
(SBX-D), and SBX-Fast.
For the germline WGS `--library-type` option, see
[Library types for germline WGS](appendix.md#library-types-for-germline-wgs).

```bash
compute_vcf_features germline-wgs \
    --library-type sbx-d \
    --genome REFERENCE.fa \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --pop-af-vcf POPAF.vcf.gz \
    --output-file VCF_FEATURES.tsv \
    --output-bed VCF_REGIONS.bed \
    --threads 16
```

For inputs, outputs, and the full CLI options, see the
[Compute VCF features](usage-compute-vcf-features.md).

## Generate a panel of normals

Aggregates BAM feature tables from normal samples into a Panel of Normals VCF.
Input: **normal-sample BAM features TSVs** → Output: **PON VCF**.

```bash
generate_pon \
    --feature-list PON_FEATURES.txt \
    --output-file pon.vcf.gz \
    --threads 16
```

For inputs, outputs, and the full CLI options, see the
[Generate a panel of normals](usage-generate-pon.md).

## Convert VCF to BED

Creates a BED file from the variant positions in a VCF, which can then be used
to restrict feature computation during model training.
Input: **VCF** → Output: **BED**.

```bash
vcf_to_bed \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --output-file vcf-regions.bed
```

For the full CLI options, see the [Convert VCF to BED](usage-vcf-to-bed.md).

## Detect sample type

Infers whether a tumor-normal Mutect2 VCF corresponds to an FFPE or a cell-line
sample. This is the standalone entry point for the same inference that
`filter_variants tumor-normal-wgs` runs under `--sample-type auto-detect`.
Input: **tumor-normal VCF** → Output: **sample type on stdout**.

```bash
detect_sample_type \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --threads 16
```

For inputs, outputs, and the full CLI options, see the
[Detect sample type](usage-detect-sample-type.md).

## Convert features to VCF

Converts a BAM features TSV into a coordinate-sorted, indexed, single-sample
VCF with synthetic genotype fields, so downstream tools that require a VCF can
read the variant sites.
Input: **BAM features TSV + reference** → Output: **single-sample VCF**.

```bash
features_to_vcf \
    --input-file BAM_FEATURES.tsv \
    --genome REFERENCE.fa \
    --sample-name SAMPLE \
    --output-file features.vcf.gz
```

For inputs, outputs, and the full CLI options, see the
[Convert features to VCF](usage-features-to-vcf.md).
