# Train models

`train_model` fits LightGBM models for variant filtering and re-genotyping from
pre-computed feature files and truth VCFs.
Compute the features first with
[compute_bam_features](usage-compute-bam-features.md) and
[compute_vcf_features](usage-compute-vcf-features.md), then train against
truth VCFs.
The resulting model files are the input to
[filter_variants](usage-filter-variants.md).

Select a use case with the subcommand: `germline-wgs`, `tumor-normal-wgs`,
or `tumor-only-wgs`.

## How it works

- ML models are trained on BAM and VCF features using [LightGBM](https://github.com/microsoft/LightGBM).
- The BAM and VCF features for each variant are converted into a numeric vector.
- Numeric vectors for all variants are appended together into a single vector, which is the input for training a LightGBM model.
- Each variant in the training data is assigned a classification label based on its genotype (or lack thereof) in the truth VCF if training germline models, or on whether the variant occurs in a positive or negative sample if somatic.
- The outputs are one or more LightGBM model file(s).

## Example command

### germline-wgs

Training has no BAM input, so pass a concrete `--library-type`.
The supported library types are SBX-Simplex (SBX-S), SBX-Duplex (SBX-D), and
SBX-Fast.
Use the corresponding `--library-type` values `sbx-s`, `sbx-d`, or `sbx-fast`.
See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs)
for the complete mapping and selection behavior.
Use `--model-type` to choose the germline WGS training defaults:
`large` is the default and uses settings for better F1 performance, but may
require more memory and run more slowly than `small`.
`small` uses settings that require less memory and may have a faster runtime.
Explicit `--snv-iterations` and `--indel-iterations` values override the
selected model-type defaults.

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

See [`train_model` overview and CLI options](#cli-options) for all available parameters.

### tumor-normal-wgs

```bash
train_model tumor-normal-wgs \
    --positive-bam-features POSITIVE_BAM_FEATURES_FILE_LIST.txt \
    --positive-vcf-features POSITIVE_VCF_FEATURES_FILE_LIST.txt \
    --negative-bam-features NEGATIVE_BAM_FEATURES_FILE_LIST.txt \
    --negative-vcf-features NEGATIVE_VCF_FEATURES_FILE_LIST.txt \
    --truth-vcfs TRUTH_VCF_FILE_LIST.txt \
    --output-file MODEL.txt.gz \
    --threads 16
```

Here the `tumor-normal-wgs` sub-command is used to train a LightGBM model from one or more samples.
Samples can be replicates with different tumor purity.

See [`train_model` overview and CLI options](#cli-options) for all available parameters.

### tumor-only-wgs

```bash
train_model tumor-only-wgs \
    --positive-bam-features POSITIVE_BAM_FEATURES_FILE_LIST.txt \
    --positive-vcf-features POSITIVE_VCF_FEATURES_FILE_LIST.txt \
    --negative-bam-features NEGATIVE_BAM_FEATURES_FILE_LIST.txt \
    --negative-vcf-features NEGATIVE_VCF_FEATURES_FILE_LIST.txt \
    --somatic-truth-vcfs SOMATIC_TRUTH_VCF_FILE_LIST.txt \
    --germline-truth-vcfs GERMLINE_TRUTH_VCF_FILE_LIST.txt \
    --output-file MODEL.txt.gz \
    --threads 16
```

The `tumor-only-wgs` sub-command trains a single-sample LightGBM model using
BAM and VCF features. Somatic and germline truth sets define positive examples,
and variants outside those truth sets define negative examples.
See [`train_model` overview and CLI options](#cli-options) for all available parameters.

## Input

Training inputs differ by use case. The `germline-wgs`, `tumor-normal-wgs`, and
`tumor-only-wgs` use cases require both BAM and VCF features, computed with
`compute_bam_features` and `compute_vcf_features`. Where multiple list files
are used, the sample order must match across them.

### germline-wgs inputs

Each training sample requires three input files, whose paths are listed one per
line in the corresponding list file:

| Input | Required | Notes |
|:---|:---|:---|
| Truth VCF | Yes | Passed with `--truth-vcfs`. High-confidence calls for the sample; the truth label per variant is read from the FORMAT `GT` field. |
| BAM features | Yes | Passed with `--positive-bam-features`. Produced by `compute_bam_features`. |
| VCF features | Yes | Passed with `--positive-vcf-features`. Produced by `compute_vcf_features`. |

The sample order must match across `--truth-vcfs`, `--positive-bam-features`,
and `--positive-vcf-features`. The command accepts one or more samples. The
selected library type determines the duplex or simplex feature schema, and
feature files from different schemas must not be mixed.
For germline whole-genome samples, high-confidence variant calls can
be obtained from the [Genome in a Bottle (GIAB) consortium](https://jimb.stanford.edu/giab/).

Example commands to compute BAM and VCF features for one training sample:

```bash
compute_vcf_features germline-wgs \
    --library-type sbx-d \
    --genome REFERENCE.fa \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --pop-af-vcf POPAF.vcf.gz \
    --threads 16 \
    --target-regions TARGET.bed \
    --output-file VCF_FEATURES.tsv \
    --output-bed VCF_REGIONS.bed

compute_bam_features germline-wgs \
    --library-type sbx-d \
    --genome REFERENCE \
    --bam-input GATK_OUTPUT \
    --threads 16 \
    --target-regions VCF_REGIONS.bed \
    --output-file BAM_FEATURES.tsv
```

The `--output-bed` argument in `compute_vcf_features` writes the variant
positions to a BED file, which is passed to `compute_bam_features` via
`--target-regions` to restrict computation to variant sites and speed up
processing. See [`compute_bam_features` overview and CLI options](usage-compute-bam-features.md#cli-options)
and [`compute_vcf_features` overview and CLI options](usage-compute-vcf-features.md#cli-options)
for more details.

List the feature TSV paths, one per line, in the order used across all list
files. The example below builds list files for three training samples:

```bash
echo /path/to/sample1/BAM_FEATURES.tsv > BAM_FEATURES_FILE_LIST.txt
echo /path/to/sample2/BAM_FEATURES.tsv >> BAM_FEATURES_FILE_LIST.txt
echo /path/to/sample3/BAM_FEATURES.tsv >> BAM_FEATURES_FILE_LIST.txt

echo /path/to/sample1/VCF_FEATURES.tsv > VCF_FEATURES_FILE_LIST.txt
echo /path/to/sample2/VCF_FEATURES.tsv >> VCF_FEATURES_FILE_LIST.txt
echo /path/to/sample3/VCF_FEATURES.tsv >> VCF_FEATURES_FILE_LIST.txt

echo /path/to/sample1/truth.vcf > TRUTH_VCF_FILE_LIST.txt
echo /path/to/sample2/truth.vcf >> TRUTH_VCF_FILE_LIST.txt
echo /path/to/sample3/truth.vcf >> TRUTH_VCF_FILE_LIST.txt
```

### tumor-normal-wgs inputs

Each training sample requires five input files, whose paths are listed one per
line in the corresponding list file:

| Input | Required | Notes |
|:---|:---|:---|
| Truth VCF | Yes | Passed with `--truth-vcfs`. High-confidence somatic calls; entries are labelled `1` (true positive). |
| Positive BAM features | Yes | Passed with `--positive-bam-features`. Positive sample BAM features. |
| Positive VCF features | Yes | Passed with `--positive-vcf-features`. Positive sample VCF features. |
| Negative BAM features | Yes | Passed with `--negative-bam-features`. Negative/healthy sample BAM features. |
| Negative VCF features | Yes | Passed with `--negative-vcf-features`. Negative/healthy sample VCF features. |

The sample order must match across `--truth-vcfs`, `--positive-bam-features`,
and `--positive-vcf-features`, and separately between `--negative-bam-features`
and `--negative-vcf-features`.

Only positive data points corresponding to a truth variant (label `1`) are used
as positives; all variants in the negative/healthy samples are labelled `0`.
Variants that appear in the truth VCF(s) are excluded from the negative training
data to prevent label conflicts. This results in an uneven distribution with
many more negatives than positives. When constructing the truth VCF, remove
low-confidence calls (for example, tumor AF < 5%) to reduce noisy labels.

Example commands to compute BAM and VCF features for one training sample:

```bash
compute_vcf_features tumor-normal-wgs \
    --genome REFERENCE.fa \
    --vcf-input GATK_OUTPUT.vcf.gz \
    --pop-af-vcf POPAF.vcf.gz \
    --output-bed VARIANT_REGIONS.bed \
    --threads 16 \
    --output-file VCF_FEATURES.tsv

compute_bam_features tumor-normal-wgs \
    --genome REFERENCE.fa \
    --bam-input GATK_OUTPUT.bam \
    --target-regions VARIANT_REGIONS.bed \
    --tumor-sample-name TUMOR_SAMPLE_NAME \
    --threads 16 \
    --output-file BAM_FEATURES.tsv
```

See [`compute_bam_features` overview and CLI options](usage-compute-bam-features.md#cli-options)
and [`compute_vcf_features` overview and CLI options](usage-compute-vcf-features.md#cli-options)
for more details.

List the feature TSV paths for positive and negative samples, one per line, in
matching order. The example below builds list files for three positive and
three negative training samples:

```bash
echo /path/to/sample1/BAM_FEATURES.tsv > POSITIVE_BAM_FEATURES_FILE_LIST.txt
echo /path/to/sample2/BAM_FEATURES.tsv >> POSITIVE_BAM_FEATURES_FILE_LIST.txt
echo /path/to/sample3/BAM_FEATURES.tsv >> POSITIVE_BAM_FEATURES_FILE_LIST.txt

echo /path/to/sample1/VCF_FEATURES.tsv > POSITIVE_VCF_FEATURES_FILE_LIST.txt
echo /path/to/sample2/VCF_FEATURES.tsv >> POSITIVE_VCF_FEATURES_FILE_LIST.txt
echo /path/to/sample3/VCF_FEATURES.tsv >> POSITIVE_VCF_FEATURES_FILE_LIST.txt

echo /path/to/sample1/truth.vcf > TRUTH_VCF_FILE_LIST.txt
echo /path/to/sample2/truth.vcf >> TRUTH_VCF_FILE_LIST.txt
echo /path/to/sample3/truth.vcf >> TRUTH_VCF_FILE_LIST.txt

echo /path/to/negative_sample1/NEGATIVE_BAM_FEATURES.tsv > NEGATIVE_BAM_FEATURES_FILE_LIST.txt
echo /path/to/negative_sample2/NEGATIVE_BAM_FEATURES.tsv >> NEGATIVE_BAM_FEATURES_FILE_LIST.txt
echo /path/to/negative_sample3/NEGATIVE_BAM_FEATURES.tsv >> NEGATIVE_BAM_FEATURES_FILE_LIST.txt

echo /path/to/negative_sample1/NEGATIVE_VCF_FEATURES.tsv > NEGATIVE_VCF_FEATURES_FILE_LIST.txt
echo /path/to/negative_sample2/NEGATIVE_VCF_FEATURES.tsv >> NEGATIVE_VCF_FEATURES_FILE_LIST.txt
echo /path/to/negative_sample3/NEGATIVE_VCF_FEATURES.tsv >> NEGATIVE_VCF_FEATURES_FILE_LIST.txt
```

### tumor-only-wgs inputs

Each training sample requires BAM features, VCF features, a somatic truth VCF,
and a germline truth VCF. Paths are listed one per line in the corresponding
list file:

| Input | Required | Notes |
|:---|:---|:---|
| Somatic truth VCF | Yes | Passed with `--somatic-truth-vcfs`. Somatic truth calls for the positive sample. |
| Germline truth VCF | Yes | Passed with `--germline-truth-vcfs`. Germline truth calls for the positive sample. |
| Positive BAM features | Yes | Passed with `--positive-bam-features`. Positive sample BAM features. |
| Positive VCF features | Yes | Passed with `--positive-vcf-features`. Positive sample VCF features. |
| Negative BAM features | No | Passed with `--negative-bam-features`. Negative sample BAM features. |
| Negative VCF features | No | Passed with `--negative-vcf-features`. Negative sample VCF features. |

The sample order must match across `--positive-bam-features`,
`--positive-vcf-features`, `--somatic-truth-vcfs`, and
`--germline-truth-vcfs`. If negative feature lists are supplied, their sample
order must also match.

## Output

The `germline-wgs` use case produces separate
SNV and indel model files; the somatic use cases (`tumor-normal-wgs` and
`tumor-only-wgs`) produce a single model file:

| Use case | Default Name(s) | Description |
|:---|:---|:---|
| germline-wgs | *snv-model.txt.gz*, *indel-model.txt.gz* | SNV model file and indel model file |
| tumor-normal-wgs, tumor-only-wgs | *model.txt.gz* | model file |

The output model file(s) must not already exist in the output directory.
If an output file name ends in `.gz`, the output is gzip-compressed; otherwise it is uncompressed.

## Tips for analysis

Expand the sections below for per-use-case training guidance.

<details>
<summary>germline-wgs tips</summary>

- The default training parameters are optimized for six samples.
  - The training parameters used in this use case were chosen to optimize training for 6 whole-genome GIAB samples, with a left-out sample being used to evaluate the model's performance.
- Adjust training parameters as needed if training using a different number of samples or number of chromosomes.

</details>

<a id="tumor-normal-wgs-tips"></a>

<details>
<summary>tumor-normal-wgs tips</summary>

- Current training parameters are arbitrarily selected and have not been optimized.
- Adjust training parameters as needed if training using a different number of samples or number of chromosomes.
- Use pairs of healthy samples (with no expected somatic variation) as negative training data.
  - Generating normal-normal pairs using the normal samples paired with a buffy coat from the tumor-normal pairs is a good option.
- Having a large imbalance between the number of negative samples relative to positive samples is expected, as all variants in a healthy normal are treated as negatives, whereas only known somatic variants are treated as positives.
- Use the `--output-training-data-tsv` option to export the assembled training feature matrix for inspection or debugging before model fitting.
- Ensure that the same settings for variant detection are applied to all training samples, whether healthy or tumor, to avoid introducing biases in the features used for model training.
  - While normal-normal pairs are not expected to have true somatic variants, they must still be processed using the same GATK Mutect2 commands to generate the BAM and VCF features for model training as tumor-normal pairs.

</details>

<a id="varied-genomic-coverages-in-training-data"></a>

<details>
<summary>Varied genomic coverages in training data</summary>

**Include a range of sequencing coverages in your training data.** The ML model learns to distinguish true somatic variants from artifacts in part from depth-related features.
If the model is trained only on samples at a single coverage, it is unlikely to generalize well to samples at higher or lower depths.
Include training samples spanning the full range of coverage you expect in your production data (e.g., 25x to 170x tumor, 25x to 127x normal for the FFPE model).
**Performance is not guaranteed for data with coverage significantly above or below the range used during training.** In particular:

- Samples with much lower coverage than the training range can show increased false negatives due to the model not having learned to handle low-depth signals.
- Samples with much higher coverage than the training range can show increased false positives or unexpected behavior, as depth-normalized features can fall outside the distribution seen during training.

Include multiple tumor purity levels in training data, as the model needs to learn to call variants across a range of tumor allele frequencies.

</details>

<details>
<summary>Cell line model training</summary>

Cell line samples are often accompanied by high-confidence, well-characterized somatic variant truth sets and are a good starting point for tumor-normal model training.

Recommended training data for a cell line model:

- **Samples**: Use well-characterized tumor cell line samples (e.g. HCC1395, H2009, HG008) paired with matched normal cell lines or blood normals.
  - Include multiple replicates to teach the model to call variants across a range of allele frequencies.
  - Use unrelated samples for training and testing to avoid overfitting.
- **Coverage**: Include training samples at varied coverage levels spanning the range expected in production data.
  - For example, the pre-trained FFPE model used between 25-170x tumor and 25-127x normal.
  - See [Varied genomic coverages in training data](#varied-genomic-coverages-in-training-data) for why coverage diversity matters.
- **Truth sets**: Use high-confidence somatic truth calls derived from orthogonal methods for cell lines with well-characterized mutational landscapes.
  - See [tumor-normal-wgs inputs](#tumor-normal-wgs-inputs) for recommended truth set quality filters.
- **Negative training data**: Use normal-normal pairs of healthy normal samples (with no expected somatic variation) as described in the general [tumor-normal-wgs tips](#tumor-normal-wgs-tips).
  - Cell line normal controls or matched blood normals are appropriate.

</details>

<details>
<summary>FFPE model training</summary>

FFPE samples introduce unique sequencing artifacts (e.g., C>T deamination artifacts, fragmentation artifacts) that require specialized training data.
Building a robust FFPE model requires careful curation of both positive and negative training samples.

Recommended approach for FFPE model training:

- **Include real FFPE tissue samples**: At least some positive training samples must come from genuine FFPE-processed tissue, not just cell line proxies.
  - FFPE artifacts are tissue- and fixation-specific, and the model must be exposed to real FFPE signal to learn to distinguish artifacts from true somatic variants.
- **Truth sets for FFPE samples**: High-confidence somatic truth calls for FFPE samples can be derived from orthogonal high-quality sequencing of the same or matched tissue (e.g., fresh-frozen tissue).
  - Apply the standard truth set quality filters described in [tumor-normal-wgs inputs](#tumor-normal-wgs-inputs).
  - In addition, be especially strict about removing low-AF calls and HP/STR variants for FFPE data, as deamination and fragmentation artifacts are more prevalent and more easily confused with real low-AF variants.
- **Include cell line samples as additional positives**: Cell lines provide reliable, deeply validated truth sets.
  - Including cell line samples as part of the positive training data strengthens the model's ability to call true somatic variants.
  - Use cell lines where the somatic truth set is well-established and high-confidence (e.g. HCC1395, H2009, HG008).
- **Paired normal requirements**: You must include **matched tumor-normal pairs where the normal sample is either blood or a non-FFPE tissue type**.
  - One of the normals in the pair must be a FFPE normal (normal tissue processed through FFPE) and another must be a blood normal from the same individual.
  - This normal-normal (NN) pairing—where one normal is FFPE and the other is blood—is used to train the model to distinguish FFPE artifacts (present in FFPE tissue but not blood) from true somatic variants.
- **Coverage range**: Include training samples across a range of coverage levels (e.g., 25x to 25x, 40x to 40x, 100x to 40x, 170x to 127x tumor-normal).
  - See [Varied genomic coverages in training data](#varied-genomic-coverages-in-training-data) for why coverage diversity matters.
- **Negative training data**: Use healthy paired normals, including FFPE normals, as negatives.
  - All variants in the negative samples are treated as artifacts or germline variation, providing the model with FFPE-specific background noise patterns to learn from.

</details>

## CLI options

Required parameters are shown in **bold**.

This command also accepts the common `-h, --help`, `-v, --version`, and
`-l, --log-level` options.
Because SVC trains different model families per workflow, the tables below keep
the use-case comparisons that affect required inputs and defaults.

### Input options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --config | Path to a configuration JSON file. The use case selected is based on the subcommand used. A file containing the default configurations for each supported use case is provided. See [Custom configurations](file-formats.md#custom-configurations) for details. | Path to an existing JSON file |

#### germline-wgs input options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--library-type** | Library type used to select germline WGS defaults. The supported library types are SBX-Simplex (SBX-S), SBX-Duplex (SBX-D), and SBX-Fast; use the corresponding values `sbx-s`, `sbx-d`, or `sbx-fast`. See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs). Required because this command has no BAM input and cannot use `auto-detect`. | One of {`sbx-s`, `sbx-d`, `sbx-fast`} |
| **--positive-bam-features** | Path to a text file that lists BAM features file path(s), one path per line. Each file in this list is generated from a single sample using the `compute_bam_features` submodule. If the sample's alignments are split into multiple BAM files, a single features file must be computed for the sample by passing all BAM files to `compute_bam_features`. The number and order of the features files must match those passed to `--positive-vcf-features` and `--truth-vcfs`. | Path to an existing text file |
| **--positive-vcf-features** | Path to a text file that lists VCF features file path(s), one path per line. Each file in this list is generated from a single sample using the `compute_vcf_features` submodule. The number and order of the features files must match those passed to `--positive-bam-features` and `--truth-vcfs`. | Path to an existing text file |
| **--truth-vcfs** | Path to a text file that lists truth VCF file path(s), one path per line. Each variant record in the VCF file must have a `GT` genotype tag. The number and order of the VCF files must match those passed to `--positive-vcf-features` and `--positive-bam-features`. | Path to an existing text file |

#### tumor-normal-wgs and tumor-only-wgs input options

These input files are shared by the somatic subcommands.
They differ in workflow-specific requirements and defaults.

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--positive-bam-features** | Path to a text file that lists BAM features file path(s), one path per line. Each file in this list is generated from a single somatic sample using the `compute_bam_features` submodule. If the sample's alignments are split into multiple BAM files, a single features file must be computed for the sample by passing all BAM files to `compute_bam_features`. | Path to an existing text file |
| --negative-bam-features | Path to a text file that lists BAM features file path(s), one path per line. Each file is generated from a single negative/healthy sample using the `compute_bam_features` submodule. | Path to an existing text file |

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--positive-vcf-features** | Path to a text file that lists VCF features file path(s), one path per line. Each file is generated from a single somatic sample using the `compute_vcf_features` submodule. The number and order must match `--positive-bam-features` and the applicable truth VCF list(s). | Path to an existing text file |
| --negative-vcf-features | Path to a text file that lists VCF features file path(s), one path per line. Each file is generated from a single negative/healthy sample using the `compute_vcf_features` submodule. The number and order must match `--negative-bam-features`. | Path to an existing text file |

#### tumor-normal-wgs input options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--truth-vcfs** | Path to a text file that lists truth VCF file path(s), one path per line. The number and order of the VCF files must match those passed to `--positive-bam-features` and `--positive-vcf-features`. | Path to an existing text file |

#### tumor-only-wgs input options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| **--somatic-truth-vcfs** | Path to a text file that lists somatic truth VCF path(s), one path per line. The number and order of the VCF files must match those passed to `--positive-bam-features`. | Path to an existing text file |
| **--germline-truth-vcfs** | Path to a text file that lists germline truth VCF path(s), one path per line. The number and order of the VCF files must match those passed to `--positive-bam-features`. | Path to an existing text file |
| --somatic-positive-ratio | Fraction of eligible somatic positive variants to keep. | Float [0.0, 1.0] inclusive; default: `1.0` |
| --germline-positive-ratio | Fraction of eligible germline positive variants to keep. | Float [0.0, 1.0] inclusive; default: `1.0` |
| --negative-ratio | Fraction of eligible negative variants to keep. | Float [0.0, 1.0] inclusive; default: `1.0` |
| --min-support | Minimum alternate-read support for a variant to be eligible for training. | Integer >= 0; default: `1` |

### Output options

#### germline-wgs output options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --snv-output-file | Output path for the SNV model. File extension `.gz` or `.zst` turns on gzip or Zstandard compression. | Non-existent file; default: `snv_model.txt.gz` |
| --indel-output-file | Output path for the indel model. File extension `.gz` or `.zst` turns on gzip or Zstandard compression. | Non-existent file; default: `indel_model.txt.gz` |
| --snv-output-training-data-tsv | Output path for the SNV model training data TSV. This file is useful for debugging purposes. | Path |
| --indel-output-training-data-tsv | Output path for the indel model training data TSV. This file is useful for debugging purposes. | Path |

#### tumor-normal-wgs and tumor-only-wgs output options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --output-file | Path to the output model file. A file name with extension `.gz` or `.zst` turns on gzip or Zstandard compression. | Non-existent file; default: `model.txt.gz` |
| --output-training-data-tsv | Path to the output model training data TSV. This file is useful for debugging purposes. | Path |

### Training options

#### germline-wgs training options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --normalize-features | Normalize feature values based on the VCF feature `normal_dp`. `filter_variants` must also use `--normalize-features` for models trained with this parameter. Option `none`: Do not normalize feature values. Option `median-dp`: Normalize read-depth related feature values by the chromosome's median DP value. | One of {`none`, `median-dp`}; default: `median-dp` |
| --model-type | Germline WGS model family whose training defaults should be used. `large` uses settings for better F1 performance, but may require more memory and run more slowly than `small`. `small` uses settings that require less memory and may have a faster runtime. | One of {`large`, `small`}; default: `large` |
| --snv-iterations | The maximum number of training iterations for the germline SNV model before stopping. Overrides the selected `--model-type` default. | Integer > 0; default: `1000` for `large`, `1500` for `small` |
| --indel-iterations | The maximum number of training iterations for the germline indel model before stopping. Overrides the selected `--model-type` default. | Integer > 0; default: `8000` for `large`, `1500` for `small` |

#### tumor-normal-wgs and tumor-only-wgs training options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --normalize-features | Normalize feature values. `filter_variants` must also use `--normalize-features` for models trained with this parameter. Option `none`: Do not normalize feature values. Option `median-dp`: Normalize read-depth related feature values by the chromosome's median DP value. | One of {`none`, `median-dp`}; default: `none` |
| --iterations | The maximum number of training iterations for the model before stopping. | Integer > 0; default: tumor-normal-wgs/tumor-only-wgs: `6000` |

### Performance options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --threads | Number of threads to use (`0` = available hardware threads). | Integer >= 0 [default: `1`] |

### Developer options

| Parameter | Description | Value(s) |
|:---|:---|:---|
| --warn-as-error | Treat warning messages as errors. | [default: `false`] |
