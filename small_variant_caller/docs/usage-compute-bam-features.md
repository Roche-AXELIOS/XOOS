# Compute BAM features

`compute_bam_features` extracts per-variant machine learning features from one
or more aligned BAM files.
It is a training-time helper: run it (together with
[compute_vcf_features](usage-compute-vcf-features.md)) to produce the feature
tables that [train_model](usage-train-model.md) fits models to.
`filter_variants` computes the same features internally at inference time, so
you do not need to run this command to filter a sample.

Select a use case with the subcommand: `germline-wgs`,
`germline-mitochondria`, `tumor-normal-wgs`, `tumor-only-wgs`,
or `panel-of-normals`.

## How it works

- One or more coordinate-sorted and indexed BAM files are read via [htslib](https://github.com/samtools/htslib) to compute ML features.
- At each alignment position, ML features are computed for both reference allele and alternate allele(s).
- Reference allele features are computed from alignment sequence matches.
- Alternate allele features are computed from mismatches, insertions, and deletions.
- Many ML features, such as the mean base quality (`baseq_mean`) of variants, can be computed.
- For a list of available BAM features, refer to the section [Model feature definition](file-formats.md#model-feature-definition)
- Feature computation can be parallelized and/or restricted to target regions via the input BED file.
- Parallelization uses [Taskflow](https://github.com/taskflow/taskflow), splitting work into per-region tasks that process only overlapping BAM records.
- The output is a tabular text file, where the rows are variants and the columns are computed features.

For the full list of BAM-derived features, see
[Model feature definition](file-formats.md#model-feature-definition).

## Example command

### germline-wgs

For this BAM-based command, `--library-type auto-detect` is the default.
The supported library types are SBX-Simplex (SBX-S), SBX-Duplex (SBX-D), and
SBX-Fast.
Use the corresponding `--library-type` values `sbx-s`, `sbx-d`, or `sbx-fast`.
See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs)
for the complete mapping and selection behavior.

```bash
compute_bam_features germline-wgs \
    --genome REFERENCE \
    --bam-input GATK_OUTPUT \
    --threads 16 \
    --target-regions VCF_REGIONS.bed \
    --output-file BAM_FEATURES.tsv
```

The `--target-regions` BED restricts feature computation to the regions of
interest; it is typically the `--output-bed` produced by
[compute_vcf_features](usage-compute-vcf-features.md).

### germline-mitochondria

```bash
compute_bam_features germline-mitochondria \
    --genome REFERENCE \
    --bam-input GATK_OUTPUT \
    --threads 16 \
    --target-regions VCF_REGIONS.bed \
    --output-file BAM_FEATURES.tsv
```

### tumor-normal-wgs

```bash
compute_bam_features tumor-normal-wgs \
    --genome REFERENCE.fa \
    --bam-input GATK_OUTPUT.bam \
    --target-regions VARIANT_REGIONS.bed \
    --tumor-sample-name TUMOR_SAMPLE_NAME \
    --threads 16 \
    --output-file BAM_FEATURES.tsv
```

### tumor-only-wgs

```bash
compute_bam_features tumor-only-wgs \
    --genome REFERENCE \
    --bam-input GATK_OUTPUT \
    --target-regions VARIANT_REGIONS.bed \
    --threads 16 \
    --output-file BAM_FEATURES.tsv
```

### panel-of-normals

```bash
compute_bam_features panel-of-normals \
    --genome REFERENCE \
    --bam-input NORMAL_SAMPLE.bam \
    --target-regions VARIANT_REGIONS.bed \
    --threads 16 \
    --output-file BAM_FEATURES.tsv
```

## Input

| Input | Required | Notes |
|:---|:---|:---|
| BAM files | Yes | Passed with `--bam-input`. One or more coordinate-sorted, indexed BAM files. |
| Reference genome | Yes | Passed with `--genome`. Reference genome FASTA and its index. |
| Target regions | No | Passed with `--target-regions`. Optional BED to restrict computation. |

To generate the required candidate BAM, see
[Generate candidate calls](generating-candidate-calls.md).

## Output

| Output | Description |
|:---|:---|
| BAM features TSV | Written to `--output-file`. Tab-separated feature table with one row per variant and one column per feature. Consumed by [train_model](usage-train-model.md). For the column definitions, see [Model feature definition](file-formats.md#model-feature-definition). |

## CLI options

Required parameters are shown in **bold**.

This command also accepts the common `-h, --help`, `-v, --version`, and
`-l, --log-level` options.
Because SVC supports multiple workflows with different defaults, the tables
below keep cross-use-case comparisons where they help the reader choose the
right parameters.

### Input options

These input options apply across workflows in addition to the use-case-specific
BAM inputs listed below.

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--genome** | Path to a reference genome (indexed) FASTA file. Reference must match the reference used by GATK HaplotypeCaller or Mutect2. | Path to an existing FASTA file (can be gzipped) |
| --config | Path to a configuration JSON file. The use case selected is based on the subcommand used. A file containing the default configurations for each supported use case is provided. See [Custom configurations](file-formats.md#custom-configurations) for details. | Path to an existing JSON file |
| --library-type | Library type used to select germline WGS defaults. The supported library types are SBX-Simplex (SBX-S), SBX-Duplex (SBX-D), and SBX-Fast; use the corresponding values `sbx-s`, `sbx-d`, or `sbx-fast`. See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs). `auto-detect` inspects every input BAM; explicit values bypass detection. | One of {`auto-detect`, `sbx-s`, `sbx-d`, `sbx-fast`} [default: `auto-detect`] |
| --skip-variants-vcf | VCF file containing variants not counted by `--max-variants-per-read` and `--max-variants-per-read-normalized`. | Path to an existing indexed VCF file |

#### Required inputs per use case

| Parameter | germline-wgs | germline-mitochondria | tumor-normal-wgs | tumor-only-wgs | panel-of-normals |
|:---|:---|:---|:---|:---|:---|
| **--bam-input** | GATK HaplotypeCaller BAM(s) | GATK Mutect2 (mitochondria mode) BAM | GATK Mutect2 BAM(s) (with read groups) | GATK Mutect2 BAM(s) | Normal-sample BAM(s) |

Input BAMs must be coordinate-sorted, indexed, and contain sequence and quality
strings for each alignment. When multiple BAMs are given they must be from the
same sample, and duplicate alignments across BAMs are counted once. For
`tumor-normal-wgs`, `--tumor-sample-name` identifies the tumor read group.

### Output options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --output-file | Path to the output features TSV file. | Non-existent file [default: `bam_features.txt`] |

### Region options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --target-regions | Restricts processing to 0-based target regions. | An existing BED file (uncompressed or gzip/BGZF-compressed `.gz`) |

### Sample metadata options

#### tumor-normal-wgs sample metadata options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--tumor-sample-name** | Name of the tumor sample's read group. Used to identify which alignments come from the tumor sample. | Sample name as a String |

### Read filtering options

These options are accepted by every use case with identical semantics; only the
default value differs by use case.
Options `--decode-yc`, `--duplex-lowbq`, and `--min-base-type` apply only when
`--sequencing-protocol` is `duplex` or `duplex-simplex`.

The **Default per use case** column uses these abbreviations.
See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs)
for the SBX-Simplex (SBX-S), SBX-Duplex (SBX-D), and SBX-Fast library-type
values.
`gwd` = `germline-wgs` with an SBX-Duplex (SBX-D) or SBX-Fast library type,
`gws` = `germline-wgs` with an SBX-Simplex (SBX-S) library type, `gm` =
`germline-mitochondria`, `tnw` = `tumor-normal-wgs`, `tow` = `tumor-only-wgs`,
and `pon` = `panel-of-normals`.

| Parameter | Description | Value(s) | Default per use case |
|:---|:---|:---|:---|
| --sequencing-protocol | The sequencing protocol used to generate the input data. Option `simplex`: raw simplex sequencing data. Option `duplex`: duplex sequencing data, e.g. SBX-Duplex (SBX-D). Option `duplex-simplex`: a hybrid of duplex and simplex sequencing data, e.g. SBX-Duplex (SBX-D) and SBX-Simplex (SBX-S). Option `umi`: UMI consensus sequencing data. | One of {`simplex`, `duplex`, `duplex-simplex`, `umi`} | gws: `simplex`<br>gwd/gm/tnw/tow/pon: `duplex` |
| --decode-yc | Decode base types for duplex consensus reads using `YC` tags. Option `none`: do not find and decode YC tags; base types are inferred from base quality. Option `consensus`: decode YC tags to infer base types for consensus reads. Option `split`: decode YC tags to split each consensus read into two constituent reads (R1, R2) before feature extraction. | One of {`none`, `consensus`, `split`} | gwd/gm/tnw/tow/pon: `consensus`<br>gws: `none` |
| --duplex-lowbq | Include or exclude `duplex_lowbq` in `duplex_dp` and `duplex_af` feature calculations. | One of {`include`, `exclude`} | gwd/tnw/tow/pon: `include`<br>gws/gm: `exclude` |
| --min-base-type | Minimum base type in duplex reads for variant support. | One of {`discordant`, `simplex`, `concordant`} | gwd/gws: `simplex`<br>gm/tnw/tow/pon: `concordant` |
| --min-mapq | Minimum mapping quality. Adjusting this parameter can alter the BAM support for variant calls. Any read alignments that support the variant but have a mapping quality less than the minimum specified will be ignored and not counted as a variant supporting alignment. | Integer [0, 255] inclusive | gwd/gws/gm/tnw/tow/pon: `1` |
| --min-bq | Minimum base quality. Adjusting this parameter can alter the BAM support for variant calls. Any read alignments that support the variant but have a base quality less than the minimum specified will be ignored and not counted as a variant supporting alignment. | Integer [0, 93] inclusive | gwd/gws/gm: `0`<br>tnw/tow/pon: `23` |
| --min-dist | Minimum distance of variant from fragment end. If the variant is found to be less than this number of bases from the end of a read's alignment, the alignment will not be counted as a variant supporting alignment. | Integer >= 0 | gwd/gws/tnw/tow/pon: `2`<br>gm: `5` |
| --min-family-size | Minimum cluster size. Raw simplex reads use a fixed family size of 1. If `--sequencing-protocol duplex`, values greater than 2 are capped at 2. | Integer >= 0 | gwd/gws/gm: `1`<br>tnw/tow/pon: `2` |
| --filter-homopolymer | Skip variant adjacent to homopolymer that spans beyond the alignment's end position. | One of {`none`, `alignment-end`} | gwd/gws: `alignment-end`<br>gm/tnw/tow/pon: `none` |
| --min-homopolymer-length | Minimum length of homopolymer in reference. Only has an effect when `--filter-homopolymer alignment-end` is set. The `--min-homopolymer-length` value is used as a threshold to determine whether or not a stretch of identical bases in a read is considered a homopolymer. Only runs of identical bases greater or equal to the value of this parameter are considered homopolymers. | Integer >= 0 | gwd/gws: `4`<br>gm: `0`<br>tnw/tow/pon: `7` |
| --max-variants-per-read | Maximum number of variants allowed per read. By default, this is set to 0, meaning there is no set maximum. If the value is set to a number greater than 0, read alignments that contain more than the value number of variants are skipped. | Integer >= 0 | gwd/gws/gm/tnw/tow/pon: `0` |
| --max-variants-per-read-normalized | Maximum number of variants allowed per read, normalized by alignment length (i.e. match, mismatch, and insertion only). By default, this is set to 0, meaning there is no set maximum. If the value is set to a number greater than 0, read alignments that contain more than the value number of variants normalized to the alignment length are skipped. | Float [0.0, 1.0] inclusive | `0` |

### Performance options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --threads | Number of threads used (`0` = available hardware threads). | Integer >= 0 [default: `1`] |
| --max-region-size-per-thread | Maximum genomic region size per thread in base pairs. | Integer > 0 [default: `16384`] |

### Developer options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --warn-as-error | Treat warning messages as errors. | [default: `false`] |
