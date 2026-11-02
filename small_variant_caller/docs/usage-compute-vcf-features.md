# Compute VCF features

`compute_vcf_features` extracts per-variant machine learning features from a
single-sample GATK VCF file.
It is a training-time helper: run it (together with
[compute_bam_features](usage-compute-bam-features.md)) to produce the feature
tables that [train_model](usage-train-model.md) fits models to.
`filter_variants` computes the same features internally at inference time, so
you do not need to run this command to filter a sample.

Select a use case with the subcommand: `germline-wgs`,
`germline-mitochondria`, `tumor-normal-wgs`, or `tumor-only-wgs`.

## How it works

- A coordinate-sorted, bgzip-compressed, and indexed VCF is read via [htslib](https://github.com/samtools/htslib) to compute ML features.
- For each VCF record, ML features are computed for both reference allele and alternate allele(s).
- Many ML features, such as the allele depth (`alt_ad`) of variants, can be computed.
- For a list of available VCF features, refer to the section [Model feature definition](file-formats.md#model-feature-definition).
- Feature computation can be parallelized and/or restricted to target regions via the input BED file.
- Parallelization uses [Taskflow](https://github.com/taskflow/taskflow), splitting work into per-region tasks that process only overlapping VCF records.
- The output is a tabular text file, where the rows are variants and the columns are computed features.

For the full list of VCF-derived features, see
[Model feature definition](file-formats.md#model-feature-definition).

## Example command

### germline-wgs

This command has no BAM input, so pass a concrete `--library-type`.
The supported library types are SBX-Simplex (SBX-S), SBX-Duplex (SBX-D), and
SBX-Fast.
Use the corresponding `--library-type` values `sbx-s`, `sbx-d`, or `sbx-fast`.
See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs)
for the complete mapping and selection behavior.

```bash
compute_vcf_features germline-wgs \
    --library-type sbx-d \
    --genome REFERENCE \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --threads 16 \
    --target-regions VCF_REGIONS.bed \
    --output-file VCF_FEATURES.tsv \
    --output-bed VCF_REGIONS.bed
```

The `--output-bed` regions can be reused as the `--target-regions` input to
[compute_bam_features](usage-compute-bam-features.md), so that both commands
compute features over the same variant loci.

### germline-mitochondria

```bash
compute_vcf_features germline-mitochondria \
    --genome REFERENCE \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --threads 16 \
    --output-file VCF_FEATURES.tsv
```

### tumor-normal-wgs

```bash
compute_vcf_features tumor-normal-wgs \
    --genome REFERENCE \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --threads 16 \
    --output-file VCF_FEATURES.tsv
```

### tumor-only-wgs

```bash
compute_vcf_features tumor-only-wgs \
    --genome REFERENCE \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --threads 16 \
    --output-file VCF_FEATURES.tsv
```

## Input

| Input | Required | Notes |
|:---|:---|:---|
| VCF file | Yes | Passed with `--vcf-input`. Coordinate-sorted, bgzip-compressed, indexed VCF file. |
| Reference genome | Yes | Passed with `--genome`. Reference genome FASTA and its index. |
| Target regions | No | Passed with `--target-regions`. Optional BED to restrict computation. |
| Population allele frequency VCF | No | Passed with `--pop-af-vcf`. Used when the corresponding features are used for scoring. |
| Panel of Normals VCF | No | Passed with `--pon-vcf`. Used when the corresponding features are used for scoring. |
| Catalog VCF | No | Passed with `--catalog-vcf`. Used when the corresponding features are used for scoring. |

To generate the required candidate VCF, see
[Generate candidate calls](generating-candidate-calls.md).

## Output

| Output | Description |
|:---|:---|
| VCF features TSV | Written to `--output-file`. Tab-separated feature table with one row per variant and one column per feature. Consumed by [train_model](usage-train-model.md). For the column definitions, see [Model feature definition](file-formats.md#model-feature-definition). |
| Regions BED | Written to `--output-bed`. Optional BED file of the extracted variant regions that can be passed to [compute_bam_features](usage-compute-bam-features.md). |

## CLI options

Required parameters are shown in **bold**.

This command also accepts the common `-h, --help`, `-v, --version`, and
`-l, --log-level` options.

### Input options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--vcf-input** | Path to an input VCF file (`*.vcf.gz`) produced by GATK Mutect2 or HaplotypeCaller. The VCF must be coordinate-sorted, bgzip-compressed, and indexed. | Path to an existing indexed VCF file |
| **--genome** | Path to a reference genome (indexed) FASTA file. Reference must match the reference used by GATK HaplotypeCaller or Mutect2. | Path to an existing FASTA file (can be gzipped) |
| --config | Path to a configuration JSON file. The use case selected is based on the subcommand used. A file containing the default configurations for each supported use case is provided. See [Custom configurations](file-formats.md#custom-configurations) for details. | Path to an existing JSON file |
| **--library-type** | Library type used to select germline WGS defaults. The supported library types are SBX-Simplex (SBX-S), SBX-Duplex (SBX-D), and SBX-Fast; use the corresponding values `sbx-s`, `sbx-d`, or `sbx-fast`. See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs). Required for `germline-wgs` because this command has no BAM input and cannot use `auto-detect`. | One of {`sbx-s`, `sbx-d`, `sbx-fast`} |

#### Feature resource options

These additional input files are shared with
[filter_variants](usage-filter-variants.md) for VCF feature extraction.

| Parameter          | Description                                                                                                                                                                                                                                                                                            | Value(s)                     |
|:--------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------|
| --interest-regions | Defines 0-based regions of interest for the VCF feature `at_interest_region`.                                                                                                                                                                                                                    | An existing BED file (uncompressed or gzip/BGZF-compressed `.gz`) |
| --pop-af-vcf       | Path to a VCF containing population allele frequency (i.e. [af-only-gnomad.hg38.vcf.gz](https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz;tab=live_object)) in the *AF* INFO field. Required only if `popaf` is a scoring feature. | Path to an existing indexed VCF file |
| --pon-vcf          | Path to an indexed Panel of Normals VCF (produced by `generate_pon`). Required only if a PON-derived feature is used for scoring.                                                                                                                                                                       | Path to an existing indexed VCF file |
| --catalog-vcf      | Path to an indexed VCF of catalog variants. When provided, the binary `in_catalog` feature is set to `1` for variants matching a catalog record (chrom + position + trimmed REF/ALT) and `0` otherwise. Required only if `in_catalog` is a scoring feature.                                             | Path to an existing indexed VCF file |

The `germline-mitochondria` workflow does not use the population allele
frequency, Panel of Normals, or catalog features.

### Required inputs by workflow

| Input | germline-wgs | germline-mitochondria | tumor-normal-wgs | tumor-only-wgs |
|:---|:---|:---|:---|:---|
| **--vcf-input** | HaplotypeCaller VCF | Mutect2 mitochondria-mode VCF | Mutect2 tumor/normal VCF | Mutect2 tumor-only VCF |
| **--genome** | Indexed reference FASTA | Indexed reference FASTA | Indexed reference FASTA | Indexed reference FASTA |
| **--pon-vcf** | Not used | Not used | Optional unless required by the selected feature schema | Required for PON-derived scoring features |

The workflow-specific feature schema determines whether optional resource VCFs
such as `--pop-af-vcf` or `--catalog-vcf` are required. See
[Custom configurations](file-formats.md#custom-configurations).

### Output options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --output-file | Path to the output features TSV file. | Non-existent file [default: `vcf_features.txt`] |
| --output-bed | Path to an output BED file for extracted features. | Path to a non-existent BED file |

#### Output BED shaping options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --left-pad | Left-padding to variant start position for the output BED file. | Integer >= 0 [default: `0`] |
| --right-pad | Right-padding to variant end position for the output BED file. | Integer >= 0 [default: `0`] |
| --collapse-distance | Distance to collapse nearby variants after padding is applied for the output BED file. | Integer >= 0 [default: `0`] |

### Region options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --target-regions | Restricts processing to 0-based target regions. | An existing BED file (uncompressed or gzip/BGZF-compressed `.gz`) |

### Performance options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --threads | Number of threads used (`0` = available hardware threads). | Integer >= 0 [default: `1`] |

### Developer options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --warn-as-error | Treat warning messages as errors. | [default: `false`] |
