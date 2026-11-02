# File formats

This page is the reference for the configuration, feature, output-field, and
SHAP TSV schemas used by the Small Variant Caller.
Usage details live on the [filter_variants](usage-filter-variants.md),
[phase_variants](usage-phase-variants.md),
[train_model](usage-train-model.md),
[compute_bam_features](usage-compute-bam-features.md),
[compute_vcf_features](usage-compute-vcf-features.md), and
[generate_pon](usage-generate-pon.md) pages.

## SHAP value TSV column definition

The optional SHAP output written by `filter_variants` contains one row per
scored variant. Its columns are:

| Column | Description |
|:---|:---|
| `VARIANT` | Canonical variant identifier. |
| `GENOTYPE` | Model genotype classification for the scored variant. |
| `PRED_ML` | Model prediction probability. |
| Feature columns | One column per model scoring feature, in model order; values are the feature's SHAP contribution. |
| `BASE_VALUE` | Model baseline value used with the feature contributions. |

The feature columns and their order are model-specific. SHAP output is intended
for targeted regions because calculating it can substantially increase runtime.

## Custom configurations

Applies to: `filter_variants`, `train_model`, `compute_bam_features`,
`compute_vcf_features`.

All SVC submodules except `vcf_to_bed`, `detect_sample_type`,
`features_to_vcf`, and `generate_pon` accept a configuration JSON file
(`--config`) containing use case-specific settings for command-line arguments,
model training parameters, and model feature lists.
Keys that correspond to command-line arguments can be overridden on the command line.

The default configuration file is `/resources/profiles_config.json` in the Docker image.
Use it as a template for custom configurations.
For the `germline-wgs` workflow, `--library-type` selects the matching defaults.
The supported library types are SBX-Simplex (SBX-S), SBX-Duplex (SBX-D), and
SBX-Fast.
Their supported `--library-type` values are `sbx-s`, `sbx-d`, and `sbx-fast`.
See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs)
for the complete mapping and behavior.
The selected library type controls defaults and feature lists. In
`train_model germline-wgs`, `--model-type` selects the training iteration counts
and LightGBM parameters unless they are overridden by explicit training options.

### Configuration JSON keys

Keys for CLI parameters:

| Key                              | Parameter                          | filter_variants | train_model | compute_bam_features | compute_vcf_features |
|:----------------------------------|:------------------------------------|:-----------------|:-------------|:----------------------|:----------------------|
| min_mapq                         | --min-mapq                         | Y               | N           | Y                    | N                    |
| min_bq                           | --min-bq                           | Y               | N           | Y                    | N                    |
| min_dist                         | --min-dist                         | Y               | N           | Y                    | N                    |
| min_family_size                  | --min-family-size                  | Y               | N           | Y                    | N                    |
| filter_homopolymer               | --filter-homopolymer               | Y               | N           | Y                    | N                    |
| min_homopolymer_length           | --min-homopolymer-length           | Y               | N           | Y                    | N                    |
| sequencing_protocol              | --sequencing-protocol              | Y               | N           | Y                    | N                    |
| decode_yc                        | --decode-yc                        | Y               | N           | Y                    | N                    |
| min_base_type                    | --min-base-type                    | Y               | N           | Y                    | N                    |
| duplex_lowbq                     | --duplex-lowbq                     | Y               | N           | Y                    | N                    |
| max_variants_per_read            | --max-variants-per-read            | N               | N           | Y                    | N                    |
| max_variants_per_read_normalized | --max-variants-per-read-normalized | N               | N           | Y                    | N                    |
| iterations                       | --iterations                       | N               | Y           | N                    | N                    |
| snv_iterations                   | --snv-iterations                   | N               | Y           | N                    | N                    |
| indel_iterations                 | --indel-iterations                 | N               | Y           | N                    | N                    |
| normalize_features               | --normalize-features               | Y               | Y           | N                    | N                    |
| max_indel_size                   | --max-indel-size                   | Y               | N           | N                    | N                    |
| min_dp_ratio                     | --min-dp-ratio                     | Y               | N           | N                    | N                    |
| snv_min_ml_score                 | --snv-min-ml-score                 | Y               | N           | N                    | N                    |
| indel_min_ml_score               | --indel-min-ml-score               | Y               | N           | N                    | N                    |
| min_tumor_support                | --min-tumor-support                | Y               | N           | N                    | N                    |
| max_normal_support               | --max-normal-support               | Y               | N           | N                    | N                    |
| min_support                      | --min-support                      | Y               | Y           | N                    | N                    |
| min_tumor_af                     | --min-tumor-af                     | Y               | N           | N                    | N                    |
| min_alt_counts                   | --min-alt-counts                   | Y               | N           | N                    | N                    |
| min_af                           | --min-af                           | Y               | N           | N                    | N                    |
| min_dp                           | --min-dp                           | Y               | N           | N                    | N                    |
| min_hom_alt_af                   | --min-hom-alt-af                   | Y               | N           | N                    | N                    |
| min_phased_af                    | --min-phased-af                    | Y               | N           | N                    | N                    |
| max_phased_af                    | --max-phased-af                    | Y               | N           | N                    | N                    |
| min_weighted_counts              | --min-weighted-counts              | Y               | N           | N                    | N                    |
| hotspot_min_weighted_counts      | --hotspot-min-weighted-counts      | Y               | N           | N                    | N                    |
| min_ml_score                     | --min-ml-score                     | Y               | N           | N                    | N                    |
| hotspot_min_ml_score             | --hotspot-min-ml-score             | Y               | N           | N                    | N                    |
| phased                           | --phased                           | Y               | N           | N                    | N                    |

Keys for non-CLI parameters:

- Primarily used in model training and feature extraction

| Key                                | Description                                              |
|:------------------------------------|:----------------------------------------------------------|
| workflow                           | Workflow name used by the configuration entry            |
| n_classes                          | Number of model classes                                  |
| bam_feature_names                  | BAM feature names                                        |
| vcf_feature_names                  | VCF feature names                                        |
| snv_scoring_names                  | SNV model scoring feature names                          |
| snv_categorical_names              | SNV model categorical feature names                      |
| indel_scoring_names                | Indel model scoring feature names                        |
| indel_categorical_names            | Indel model categorical feature names                    |
| scoring_names                      | Scoring feature names                                    |
| categorical_names                  | Categorical feature names                                |
| use_vcf_features                   | Whether to use VCF features in model training            |
| model_thresholds                   | Model-specific filtering thresholds keyed by sample type and aligner |
| snv_model_file                     | SNV model file name                                      |
| indel_model_file                   | Indel model file name                                    |
| model_lgbm_params                  | LightGBM parameters for model training                   |
| model_lgbm_prediction_params       | LightGBM parameters for model prediction/inference       |
| snv_model_lgbm_params              | LightGBM parameters for SNV model training               |
| snv_model_lgbm_prediction_params   | LightGBM parameters for SNV model prediction/inference   |
| indel_model_lgbm_params            | LightGBM parameters for indel model training             |
| indel_model_lgbm_prediction_params | LightGBM parameters for indel model prediction/inference |
| germline_model_type_training_defaults | Nested `large` and `small` training defaults used by `train_model germline-wgs --model-type` |

Back to [train_model input](usage-train-model.md#input),
[compute_bam_features input](usage-compute-bam-features.md#input), and
[compute_vcf_features input](usage-compute-vcf-features.md#input).

## Added and modified fields in the output VCF

Applies to: `filter_variants`.

`filter_variants` adds or modifies fields in the output VCF as part of the
filtering and re-genotyping process. The fields depend on the use case. See
[Understanding PRED_ML, GQ, and QUAL in the output VCF](appendix.md#understanding-pred_ml-gq-and-qual-in-the-output-vcf)
for interpretation details.

### Duplex allelic-depth fields

These FORMAT fields are present only when the sequencing protocol is duplex.
They apply to duplex `germline-wgs` runs using SBX-Duplex (SBX-D) or SBX-Fast,
and to `tumor-normal-wgs` and `tumor-only-wgs` runs.
See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs)
for the `--library-type` values that select duplex processing.

| Field | Type   | Definition                                                                                  |
|:-------|:--------|:---------------------------------------------------------------------------------------------|
| ADC   | FORMAT | Allelic depths of concordant molecules for the ref and alt alleles in the order listed.     |
| ADS   | FORMAT | Allelic depths of simplex molecules for the ref and alt alleles in the order listed.        |
| ADD   | FORMAT | Allelic depths of discordant molecules for the ref and alt alleles in the order listed.     |
| ADL   | FORMAT | Allelic depths of lowbq molecules for the ref and alt alleles in the order listed.          |

### germline-wgs

The germline use cases share the same output fields.

| Field         | Type   | Definition                                                                                                                                                                                                                                                             |
|:---------------|:--------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| QUAL          | QUAL   | Variant quality score. The value ranges from 0 to 255. It is calculated as `-10 * log10(1 - P(variant))`, where `P(variant)` is the sum of the prediction probabilities of all non-reference genotypes.                                                               |
| ML_PROCESSED  | FORMAT | Indicates whether a variant was processed by `filter_variants` ML-based filtering. Assigned a value of `1` if processed, or `0` otherwise. A record can be failed due to other criteria prior to ML filtering and thus have a filter value of `FAIL` but not be ML processed. |
| PRED_ML       | FORMAT | Prediction probability score for the predicted genotype. The value ranges from 0 to 1.                                                                                                                                                                               |
| GQ            | FORMAT | Genotype quality. The value ranges from 0 to 255. It is calculated as `-10 * log10(1 - P(genotype))`, where `P(genotype)` is the prediction probability score for the variant genotype. `GQ` is essentially the Phred quality score for the value of `PRED_ML`.        |
| DP            | FORMAT | The sum of REF allele support and all ALT allele support at the variant position computed based on the input BAM file.                                                                                                                                                |
| GT            | FORMAT | The re-genotyped genotype (`GT`) value assigned to the record by the `filter_variants` submodule using the ML models.                                                                                                                                                 |
| AD            | FORMAT | The REF and ALT support for the variant computed from the input BAM file. Each comma-separated value is the sum of non-discordant supporting reads and 0.5 x discordant supporting reads, rounded to the nearest integer.                                              |
| GNOMAD_AF     | FORMAT | The gnomAD allele frequency (AF) for this variant if present within the specified `--pop-af-vcf` VCF file. If the gnomAD AF is not used by the ML filtering or is missing, a value of `-1` is assigned.                                                                 |
| REF_AVG_MAPQ  | FORMAT | The average mapping quality for reads supporting the REF allele. If not used in the ML filtering process or is missing, a value of `-1` is assigned.                                                                                                                   |
| ALT_AVG_MAPQ  | FORMAT | The average mapping quality for reads supporting an ALT allele. If not used in the ML filtering process or is missing, a value of `-1` is assigned.                                                                                                                    |
| REF_AVG_DIST  | FORMAT | The average distance from the variant site to the nearest end of the read for reads supporting the REF allele. If not used in the ML filtering process or is missing, a value of `-1` is assigned.                                                                     |
| ALT_AVG_DIST  | FORMAT | The average distance from the variant site to the nearest end of the read for reads supporting an ALT allele. If not used in the ML filtering process or is missing, a value of `-1` is assigned.                                                                      |
| DENSITY_100BP | FORMAT | The number of other variants called by GATK HaplotypeCaller found within 100bp of the variant site. If not used in the ML filtering process or is missing, a value of `-1` is assigned.                                                                                |

Duplex allelic-depth fields (`ADC`, `ADS`, `ADD`, `ADL`) are also present for duplex data. See [Duplex allelic-depth fields](#duplex-allelic-depth-fields).

### germline-mitochondria

The mitochondrial workflow filters Mutect2 candidate records using allele
frequency, depth, and feature availability. It does not apply a germline-WGS
machine-learning genotype model.

| Field | Type | Definition |
|:---|:---|:---|
| GT | FORMAT | Genotype assigned from allele frequency using `--min-hom-alt-af`. |
| AD | FORMAT | Reference and alternate support from BAM evidence. |
| AF | FORMAT | Alternate allele frequency from BAM evidence. |
| DP | FORMAT | Total depth from BAM evidence. |

The workflow can add `PASS`, `low_af`, `low_dp`, `missing_feature`, and
`blocklist` FILTER values. `QUAL`, `GQ`, and `PRED_ML` are not ML score fields
for this workflow.

### tumor-normal-wgs

| Field   | Type   | Definition                                                                                                                                                  |
|:---------|:--------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| QUAL    | QUAL   | Variant quality score. The value ranges from 0 to 255. It is calculated as `-10 * log10(1 - PRED_ML)`, where `PRED_ML` is the prediction probability score. |
| GT      | FORMAT | Genotype. Set to `0/0` for the normal sample and `0/1` for the tumor sample on passing variants.                                                            |
| AD      | FORMAT | Allele depth. One value per allele per sample (Number=R): `[ref_depth, alt_depth]` for each sample.                                                         |
| AF      | FORMAT | Allele frequency per sample, computed from BAM evidence.                                                                                                    |
| DP      | FORMAT | Total read depth per sample.                                                                                                                                |
| GQ      | FORMAT | Genotype quality per sample.                                                                                                                                |
| PRED_ML | FORMAT | Prediction probability score for the variant. The value ranges from 0 to 1.                                                                                 |
| BQ      | FORMAT | The mean base quality of the alt allele.                                                                                                                    |
| MQ      | FORMAT | The mean mapping quality of the alt allele.                                                                                                                 |
| DT      | FORMAT | The mean distance from the end of the alignment of the alt allele.                                                                                          |
| REFBQ   | INFO   | Mean base quality of the reference allele.                                                                                                                  |
| REFMQ   | INFO   | Mean mapping quality of the reference allele.                                                                                                               |
| ALTBQ   | INFO   | Mean base quality of the alt allele.                                                                                                                        |
| ALTMQ   | INFO   | Mean mapping quality of the alt allele.                                                                                                                     |
| SUBTYPE | INFO   | Numeric value representing the type of substitution or if variant is an indel.                                                                              |
| CONTEXT | INFO   | Base context around the variant, one upstream and one downstream base.                                                                                      |
| POPAF   | INFO   | The population allele frequency for this variant, 0 if variant not found in provided pop-af-vcf or no pop-af-vcf file provided.                             |
| NALOD   | INFO   | Normal artifact log odds score (passed through from GATK Mutect2).                                                                                          |
| NLOD    | INFO   | Normal log odds score (passed through from GATK Mutect2).                                                                                                   |
| TLOD    | INFO   | Tumor log odds score (passed through from GATK Mutect2).                                                                                                    |
| MPOS    | INFO   | Median position of the variant in supporting reads (passed through from GATK Mutect2).                                                                      |

Duplex allelic-depth fields (`ADC`, `ADS`, `ADD`, `ADL`) are also present for duplex data. See [Duplex allelic-depth fields](#duplex-allelic-depth-fields).

#### FILTER tags in the tumor-normal-wgs output VCF

The `tumor-normal-wgs` output VCF contains the following FILTER tags to describe the filtering status of each variant record.

| FILTER tag           | Description [default]                                            |
|:----------------------|:------------------------------------------------------------------|
| `PASS`               | Site contains at least one allele that passes filters            |
| `ML_Score`           | Filtered due to low ML model prediction score                    |
| `max_indel_size`     | Filtered due to indel size exceeding maximum allowed size [`30`] |
| `max_normal_support` | Filtered due to high supporting reads in normal sample [`2`]     |
| `min_dp_ratio`       | Filtered due to low DP relative to chromosome median DP [`0.15`] |
| `min_tumor_af`       | Filtered due to low allele frequency in tumor sample [`0.05`]    |
| `min_tumor_support`  | Filtered due to low supporting reads in tumor sample [`1`]       |
| `missing_feature`    | Filtered due to missing ML feature                               |
| `non_acgt_ref_alt`   | Filtered due to non-ACGT reference or alt allele                 |

### tumor-only-wgs

The `tumor-only-wgs` output is a single-sample somatic VCF. The model score is
stored in `PRED_ML`, and passing records receive a fixed `GQ` of 100 because the
binary model does not predict genotype likelihoods.

| Field   | Type   | Definition |
|:---------|:--------|:------------|
| QUAL    | QUAL   | Variant quality calculated as `-10 * log10(1 - PRED_ML)` when ML scoring runs. Early failures (`non_acgt_ref_alt`, `max_indel_size`, and `missing_feature`) have `QUAL=0`; hard-filter failures retain the default score, and `ML_Score` failures retain the model-derived score. |
| AD      | FORMAT | Reference and alternate support from BAM evidence. |
| AF      | FORMAT | Alternate allele frequency from BAM evidence. |
| DP      | FORMAT | Total reference plus alternate support from BAM evidence. |
| GQ      | FORMAT | Fixed at `100` for records that reach hard filtering or ML scoring, including `min_support`, `min_af`, and `ML_Score` failures. Early failures (`non_acgt_ref_alt`, `max_indel_size`, and `missing_feature`) have `GQ=0`. |
| PRED_ML | FORMAT | Binary model probability for the candidate variant. |
| BQ      | FORMAT | Mean base quality of the alternate allele. |
| MQ      | FORMAT | Mean mapping quality of the alternate allele. |
| DT      | FORMAT | Mean distance from the alignment end for the alternate allele. |
| REFBQ   | INFO   | Mean base quality of the reference allele. |
| REFMQ   | INFO   | Mean mapping quality of the reference allele. |
| ALTBQ   | INFO   | Mean base quality of the alternate allele. |
| ALTMQ   | INFO   | Mean mapping quality of the alternate allele. |
| TLOD    | INFO   | Tumor log odds score passed through from GATK Mutect2. |
| MPOS    | INFO   | Median position of the variant in supporting reads. |
| SUBTYPE | INFO   | Numeric substitution or indel subtype. |
| CONTEXT | INFO   | Sequence context around the variant. |
| POPAF   | INFO   | Population allele frequency used by the feature schema. |

For duplex `tumor-only-wgs` runs, `ADC`, `ADS`, `ADD`, and `ADL` are also
present as `Number=R` FORMAT fields. See [Duplex allelic-depth fields](#duplex-allelic-depth-fields).

#### FILTER tags in the tumor-only-wgs output VCF

The `tumor-only-wgs` output VCF contains the following FILTER tags:

| FILTER tag | Description |
|:---|:---|
| `PASS` | Site contains at least one allele that passes filters. |
| `ML_Score` | Filtered due to a low ML model prediction score. |
| `max_indel_size` | Filtered due to an indel size greater than 30 bp. |
| `min_af` | Filtered due to allele frequency below the model-specific threshold. |
| `min_support` | Filtered due to low alternative allele support. |
| `missing_feature` | Filtered due to a missing ML feature. |
| `non_acgt_ref_alt` | Filtered due to a non-ACGT reference or alternative allele. |

Back to [filter_variants output](usage-filter-variants.md#output).

## Variant phaser VCF fields

Applies to: `phase_variants`.

`phase_variants` adds phasing fields to records in the output VCF.
When merging is enabled, it also adds an INFO field and can add a FILTER tag
for the merged representation and its component records.

| Field | Type | Definition |
|:---|:---|:---|
| PGT | FORMAT | Physical phasing genotype for the tumor read group. `0\|1` and `1\|0` identify the two haplotype assignments. |
| PID | FORMAT | Physical phasing identifier shared by records in the same phase block. |
| PS | FORMAT | Phase-set identifier, represented by the 1-based position of the leftmost variant in the block. |
| MNVTAG | INFO | Comma-separated list of the component variants represented by a phased block. The field is emitted when `--report-mnv` is set and is enabled automatically for records affected by merging. |
| PHASED_COMPONENT | FILTER | Original component record of a merged phase block when `--keep-components` is enabled. The merged record is emitted separately with `PASS`. |

Back to [phase_variants output](usage-phase-variants.md#output).

## Panel of Normals VCF

Applies to: `generate_pon`, `compute_vcf_features`, `filter_variants`.

`generate_pon` writes a VCF containing variants observed across the supplied
normal-sample BAM feature tables.
Each record has `FILTER=PASS`, no sample FORMAT data, and the following INFO
fields:

| Field | Type | Description |
|:---|:---|:---|
| `PON_OCC` | Integer | Number of input samples containing the variant. |
| `PON_OCC_RATIO` | Float | Fraction of input samples containing the variant. |
| `PON_MEAN_DUPLEX_AF` | Float | Mean `duplex_af` across samples containing the variant. |
| `PON_MIN_DUPLEX_AF` | Float | Minimum `duplex_af` across samples containing the variant. |
| `PON_MAX_DUPLEX_AF` | Float | Maximum `duplex_af` across samples containing the variant. |

Back to [generate_pon output](usage-generate-pon.md#output),
[compute_vcf_features input](usage-compute-vcf-features.md#input), and
[filter_variants input](usage-filter-variants.md#input).

## QC metrics output file

Applies to: `filter_variants`.

Every `filter_variants` workflow writes a tab-separated QC metrics summary to `metrics/metrics.tsv` inside the output
directory. Use `--output-metrics` to change the path (relative to `--output-dir`). The file uses a key-value layout with
columns `metric_name`, `value`, `percentage`, and `denominator`, preceded by a `##RocheCommandLine` metadata line and a
`##Sample` line. A `percentage` or ratio that cannot be computed (for example, when the denominator is zero) is reported
as `NA`.

Counts are expressed in terms of variants, meaning ALT alleles for which both REF and ALT contain only A, C, G, or T, rather than VCF records.
A multi-allelic record with two such ALT alleles contributes two to the variant totals.
The genotype breakdown is also variant-level, so a `1/2` site with two ACGT-only ALT alleles contributes two multi-allelic heterozygous variants.

| Metric                | Description |
|:-----------------------|:-------------|
| qc_status             | `pass` when there is at least one output PASS variant, otherwise `fail`. |
| predicted_sex         | Biological sex predicted: `male`, `female`, or `unknown`. |
| sample_type           | Sample type specified or detected (`ffpe` or `cell-line`); `NA` if not applicable. |
| num_input_variants    | Number of variants in the input VCF. |
| num_output_variants   | Number of variants in the output VCF. |
| num_pass_variants     | Number of output variants on records with `FILTER=PASS`. |
| num_non_pass_variants | Number of output variants on records not marked `PASS`. |
| num_pass_snvs         | Number of PASS SNVs. |
| num_pass_insertions   | Number of PASS insertions. |
| num_pass_deletions    | Number of PASS deletions. |
| num_pass_biallelic    | PASS variants on records with a bi-allelic heterozygous genotype (`GT=0/1`). |
| num_pass_multiallelic | PASS variants on records with a multi-allelic heterozygous genotype (`GT=1/2`). |
| num_pass_homozygous   | PASS variants on records with a homozygous ALT genotype (`GT=1/1` or haploid `GT=1`). |
| num_pass_heterozygous | PASS variants on records with any heterozygous genotype (`GT=0/1` or `GT=1/2`). |
| het_hom_ratio         | Ratio of heterozygous to homozygous PASS variants. |
| num_transitions       | Number of PASS SNVs that are transitions (A↔G or C↔T). |
| num_transversions     | Number of PASS SNVs that are transversions (a purine↔pyrimidine change: A↔C, A↔T, G↔C, or G↔T). |
| ti_tv_ratio           | Ratio of transitions to transversions. |
| `num_filter_<TAG>`    | Number of non-PASS variants carrying FILTER tag `<TAG>`. One row per distinct tag. |

The metrics satisfy the following relationships:

- `num_input_variants = num_output_variants`
- `num_output_variants = num_pass_variants + num_non_pass_variants`
- `num_pass_variants = num_pass_snvs + num_pass_insertions + num_pass_deletions`
- `num_pass_variants = num_pass_heterozygous + num_pass_homozygous`
- `num_pass_heterozygous = num_pass_biallelic + num_pass_multiallelic`
- `num_pass_snvs = num_transitions + num_transversions`
- `het_hom_ratio = num_pass_heterozygous / num_pass_homozygous`
- `ti_tv_ratio = num_transitions / num_transversions`

The het/hom and Ti/Tv ratios are most meaningful for germline workflows; they are still reported for other workflows for
informational purposes.

Back to [filter_variants output](usage-filter-variants.md#output).

## Model feature definition

Applies to: `compute_bam_features`, `compute_vcf_features`.

This section defines the tabular feature schemas produced during training-time
feature extraction.

### Variant ID features

Features for identifying variants, used by both BAM and VCF features.

| Name          | Description                                                                                                                                                     |
|:--------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| chrom         | Chromosome name                                                                                                                                                 |
| pos           | 1-based start position of variant in the chromosome                                                                                                             |
| ref           | Reference (REF) allele                                                                                                                                          |
| alt           | Alternate (ALT) allele                                                                                                                                          |
| alt_len       | Length of ALT-allele minus length of REF-allele                                                                                                                 |
| subtype_index | Substitution type index; `1: A>C`, `2: A>G`, `3: A>T`, `4: C>A`, `5: C>G`, `6: C>T`, `7: G>A`, `8: G>C`, `9: G>T`, `10: indel`, `11: T>A`, `12: T>C`, `13: T>G` |
| variant_type  | `0`:snv, `1`:deletion, `2`:insertion, `3`:unknown                                                                                                               |

### BAM features

Features for both ALT- and REF-alleles extracted from alignment records in
BAM files.

NOTE:

- "lowbq" features (with `_lowbq` suffix) are designated for reads with
  non-ACGT characters or low base quality in the homopolymer base(s) at or
  nearby the variant's start position.
  - The `--duplex-lowbq` option (or `duplex_lowbq` in config) controls whether
    "lowbq" duplex reads are included in the calculation of `duplex_dp` and
    `duplex_af` features.
- "nonhp" features (with `ref_nonhp_` prefix) are designated for reads
  supporting the REF-allele, where the sequence from the REF-allele's start
  position to the alignment end position is not a homopolymer.
- Any BAM feature can be computed for tumor and normal samples by attaching the
  `tumor_` and `normal_` prefixes, respectively.
  - Example: `tumor_duplex_af` and `normal_duplex_af` are derived from
    `duplex_af` for tumor and normal samples, respectively.
  - `bam_tn_af_ratio` is computed using `tumor_duplex_af` and
    `normal_duplex_af`.
    - So, `tumor_bam_tn_ratio` and `normal_bam_tn_ratio` have the same value as
      `bam_tn_af_ratio`.

| Name                      | Description                                                                                                                              |
|:--------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------|
| adt                       | Absolute value of (`duplex` - `ref_nonhp_support`) / (`duplex` + `ref_nonhp_support`);                                                   |
| adtl                      | Log 10 of (`duplex` + `ref_nonhp_support`) \* `adt`                                                                                      |
| alignmentbias             | Alignment bias; 0.25 - `f` \* (1 - `f`), where `f` is the percentage of reads aligned to the reverse strand of the reference genome      |
| bam_tn_af_ratio           | Ratio of tumor sample AF to the sum of tumor and normal sample AF; `tumor_duplex_af` / (`tumor_duplex_af` + `normal_duplex_af`)          |
| baseq_af                  | Allele frequency based on base quality for the ALT-allele                                                                                |
| baseq_lt20_count          | Number of reads supporting ALT-allele with base quality less than 20                                                                     |
| baseq_lt20_ratio          | Proportion of reads supporting ALT-allele with base quality less than 20                                                                 |
| baseq_max                 | Maximum base quality in read base(s) supporting the ALT-allele                                                                           |
| baseq_mean                | Mean base quality in read base(s) supporting the ALT-allele                                                                              |
| baseq_min                 | Minimum base quality in read base(s) supporting the ALT-allele                                                                           |
| baseq_sum                 | Sum of base quality in read base(s) supporting the ALT-allele                                                                            |
| context                   | Two characters representing the base before and after the variant's start position                                                       |
| context_index             | `context` represented as an integer                                                                                                      |
| distance_max              | Maximum distance from the variant's start position to the nearest alignment end                                                          |
| distance_mean             | Mean distance from the variant's start position to the nearest alignment end                                                             |
| distance_mean_lowbq       | Mean distance from the variant's start position to the nearest alignment end in "lowbq" reads only                                       |
| distance_min              | Minimum distance from the variant's start position to the nearest alignment end                                                          |
| distance_sum              | Sum of distances from the variant's start position to the nearest alignment end                                                          |
| distance_sum_lowbq        | Sum of distances from the variant's start position to the nearest alignment end in "lowbq" reads only                                    |
| duplex                    | Number of duplex reads supporting the ALT-allele                                                                                         |
| duplex_af                 | Allele frequency based on duplex reads supporting the ALT-allele                                                                         |
| duplex_dp                 | Depth of duplex reads supporting the ALT-allele                                                                                          |
| duplex_lowbq              | Number of duplex reads supporting the ALT-allele in "lowbq" reads only                                                                   |
| familysize_lt3_count      | Number of reads supporting ALT-allele with family size less than 3                                                                       |
| familysize_lt3_ratio      | Proportion of reads supporting ALT-allele with family size less than 3                                                                   |
| familysize_lt5_count      | Number of reads supporting ALT-allele with family size less than 5                                                                       |
| familysize_lt5_ratio      | Proportion of reads supporting ALT-allele with family size less than 5                                                                   |
| familysize_mean           | Mean family size of reads supporting the ALT-allele                                                                                      |
| familysize_sum            | Sum of family size of reads supporting the ALT-allele                                                                                    |
| indel_af                  | `duplex` / (`duplex` + `ref_nonhp_support`)                                                                                              |
| mapq_af                   | Allele frequency based on mapping quality                                                                                                |
| mapq_lt20_count           | Number of reads supporting ALT-allele with mapping quality less than 20                                                                  |
| mapq_lt20_ratio           | Proportion of reads supporting ALT-allele with mapping quality less than 20                                                              |
| mapq_lt30_count           | Number of reads supporting ALT-allele with mapping quality less than 30                                                                  |
| mapq_lt30_ratio           | Proportion of reads supporting ALT-allele with mapping quality less than 30                                                              |
| mapq_lt40_count           | Number of reads supporting ALT-allele with mapping quality less than 40                                                                  |
| mapq_lt40_ratio           | Proportion of reads supporting ALT-allele with mapping quality less than 40                                                              |
| mapq_lt60_count           | Number of reads supporting ALT-allele with mapping quality less than 60                                                                  |
| mapq_lt60_ratio           | Proportion of reads supporting ALT-allele with mapping quality less than 60                                                              |
| mapq_max                  | Maximum mapping quality in read base(s) supporting the ALT-allele                                                                        |
| mapq_mean                 | Mean mapping quality in read base(s) supporting the ALT-allele                                                                           |
| mapq_mean_lowbq           | Mean mapping quality in read base(s) supporting the ALT-allele in "lowbq" reads only                                                     |
| mapq_min                  | Minimum mapping quality in read base(s) supporting the ALT-allele                                                                        |
| mapq_sum                  | Sum of mapping quality in read base(s) supporting the ALT-allele                                                                         |
| mapq_sum_lowbq            | Sum of mapping quality in read base(s) supporting the ALT-allele in "lowbq" reads only                                                   |
| minusonly                 | Number of reads supporting the ALT-allele on the minus strand                                                                            |
| nonduplex                 | Number of non-duplex reads supporting the ALT-allele                                                                                     |
| normal_baseq_mean         | Mean base quality of reads supporting the ALT-allele in normal sample only                                                               |
| normal_distance_mean      | Mean distance from the variant's start position to the nearest alignment end of reads supporting the ALT-allele in normal sample only    |
| normal_mapq_mean          | Mean mapping quality of reads supporting the ALT-allele in normal sample only                                                            |
| normal_ref_support        | Number of reads supporting the REF-allele in normal sample only                                                                          |
| normal_support            | Number of reads supporting the ALT-allele in normal sample only                                                                          |
| num_alt                   | Number of ALT-alleles at the corresponding REF-allele's start position                                                                   |
| plusonly                  | Number of reads supporting the ALT-allele on the plus strand                                                                             |
| ref_baseq_af              | Allele frequency based on base quality for the REF-allele                                                                                |
| ref_baseq_lt20_count      | Number of reads supporting REF-allele with base quality less than 20                                                                     |
| ref_baseq_lt20_ratio      | Proportion of reads supporting REF-allele with base quality less than 20                                                                 |
| ref_baseq_mean            | Mean base quality in read base(s) supporting the REF-allele                                                                              |
| ref_baseq_sum             | Sum of base quality in read base(s) supporting the REF-allele                                                                            |
| ref_distance_max          | Maximum distance from the variant's start position to the nearest alignment end of reads supporting the REF-allele                       |
| ref_distance_mean         | Mean distance from the variant's start position to the nearest alignment end of reads supporting the REF-allele                          |
| ref_distance_mean_lowbq   | Mean distance from the variant's start position to the nearest alignment end of reads supporting the REF-allele in "lowbq" reads only    |
| ref_distance_min          | Minimum distance from the variant's start position to the nearest alignment end of reads supporting the REF-allele                       |
| ref_distance_sum          | Sum of distances from the variant's start position to the nearest alignment end of reads supporting the REF-allele                       |
| ref_distance_sum_lowbq    | Sum of distances from the variant's start position to the nearest alignment end of reads supporting the REF-allele in "lowbq" reads only |
| ref_duplex_af             | Allele frequency based on duplex reads supporting the REF-allele                                                                         |
| ref_duplex_lowbq          | Number of duplex reads supporting the REF-allele in "lowbq" reads only                                                                   |
| ref_familysize_lt3_count  | Number of reads supporting REF-allele with family size less than 3                                                                       |
| ref_familysize_lt3_ratio  | Proportion of reads supporting REF-allele with family size less than 3                                                                   |
| ref_familysize_lt5_count  | Number of reads supporting REF-allele with family size less than 5                                                                       |
| ref_familysize_lt5_ratio  | Proportion of reads supporting REF-allele with family size less than 5                                                                   |
| ref_familysize_mean       | Mean family size of reads supporting the REF-allele                                                                                      |
| ref_familysize_sum        | Sum of family size of reads supporting the REF-allele                                                                                    |
| ref_mapq_af               | Allele frequency based on mapping quality for the REF-allele                                                                             |
| ref_mapq_lt20_count       | Number of reads supporting REF-allele with mapping quality less than 20                                                                  |
| ref_mapq_lt20_ratio       | Proportion of reads supporting REF-allele with mapping quality less than 20                                                              |
| ref_mapq_lt30_count       | Number of reads supporting REF-allele with mapping quality less than 30                                                                  |
| ref_mapq_lt30_ratio       | Proportion of reads supporting REF-allele with mapping quality less than 30                                                              |
| ref_mapq_lt40_count       | Number of reads supporting REF-allele with mapping quality less than 40                                                                  |
| ref_mapq_lt40_ratio       | Proportion of reads supporting REF-allele with mapping quality less than 40                                                              |
| ref_mapq_lt60_count       | Number of reads supporting REF-allele with mapping quality less than 60                                                                  |
| ref_mapq_lt60_ratio       | Proportion of reads supporting REF-allele with mapping quality less than 60                                                              |
| ref_mapq_max              | Maximum mapping quality in read base(s) supporting the REF-allele                                                                        |
| ref_mapq_mean             | Mean mapping quality in read base(s) supporting the REF-allele                                                                           |
| ref_mapq_mean_lowbq       | Mean mapping quality in read base(s) supporting the REF-allele in "lowbq" reads only                                                     |
| ref_mapq_min              | Minimum mapping quality in read base(s) supporting the REF-allele                                                                        |
| ref_mapq_sum              | Sum of mapping quality in read base(s) supporting the REF-allele                                                                         |
| ref_mapq_sum_lowbq        | Sum of mapping quality in read base(s) supporting the REF-allele in "lowbq" reads only                                                   |
| ref_nonhp_baseq_max       | Maximum base quality in read base(s) supporting the REF-allele in "nonhp" reads only                                                     |
| ref_nonhp_baseq_mean      | Mean base quality in read base(s) supporting the REF-allele in "nonhp" reads only                                                        |
| ref_nonhp_baseq_min       | Minimum base quality in read base(s) supporting the REF-allele in "nonhp" reads only                                                     |
| ref_nonhp_baseq_sum       | Sum of base quality in read base(s) supporting the REF-allele in "nonhp" reads only                                                      |
| ref_nonhp_mapq_lt20_count | Number of reads supporting REF-allele with mapping quality less than 20 in "nonhp" reads only                                            |
| ref_nonhp_mapq_lt20_ratio | Proportion of reads supporting REF-allele with mapping quality less than 20 in "nonhp" reads only                                        |
| ref_nonhp_mapq_lt30_count | Number of reads supporting REF-allele with mapping quality less than 30 in "nonhp" reads only                                            |
| ref_nonhp_mapq_lt30_ratio | Proportion of reads supporting REF-allele with mapping quality less than 30 in "nonhp" reads only                                        |
| ref_nonhp_mapq_lt40_count | Number of reads supporting REF-allele with mapping quality less than 40 in "nonhp" reads only                                            |
| ref_nonhp_mapq_lt40_ratio | Proportion of reads supporting REF-allele with mapping quality less than 40 in "nonhp" reads only                                        |
| ref_nonhp_mapq_lt60_count | Number of reads supporting REF-allele with mapping quality less than 60 in "nonhp" reads only                                            |
| ref_nonhp_mapq_lt60_ratio | Proportion of reads supporting REF-allele with mapping quality less than 60 in "nonhp" reads only                                        |
| ref_nonhp_mapq_max        | Maximum mapping quality in read base(s) supporting the REF-allele in "nonhp" reads only                                                  |
| ref_nonhp_mapq_mean       | Mean mapping quality in read base(s) supporting the REF-allele in "nonhp" reads only                                                     |
| ref_nonhp_mapq_min        | Minimum mapping quality in read base(s) supporting the REF-allele in "nonhp" reads only                                                  |
| ref_nonhp_mapq_sum        | Sum of mapping quality in read base(s) supporting the REF-allele in "nonhp" reads only                                                   |
| ref_nonhp_support         | Number of reads supporting the REF-allele in "nonhp" reads only                                                                          |
| ref_nonhp_weighted_depth  | Weighted depth of reads supporting the REF-allele in "nonhp" reads only; sum of ((`baseq` / 138) \* (`mapq` / 60)) in all "nonhp" reads  |
| ref_support               | Number of reads supporting the REF-allele                                                                                                |
| ref_weighted_depth        | Weighted depth of reads supporting the REF-allele; sum of ((`baseq` / 138) * (`mapq` / 60)) in all reads supporting the REF-allele       |
| strandbias                | Strand bias; `f` \* (1 - `f`) where `f` = (`duplex` \* 0.5 + `plusonly`) / (`duplex` + `plusonly` + `minusonly`)                         |
| support                   | Number of reads supporting the ALT-allele                                                                                                |
| support_reverse           | Number of reads supporting the ALT-allele on the reverse strand                                                                          |
| weighted_depth            | Weighted depth of reads supporting the ALT-allele; sum of ((`baseq` / 138) \* (`mapq` / 60)) in all reads supporting the ALT-allele      |
| weightedscore             | Weighed score of duplex and non-duplex counts                                                                                            |

### VCF features

Features for both ALT- and REF-alleles extracted from the VCF file, reference
genome, and other associated resources, such as population allele frequencies.

NOTE:

- Features for short tandem repeats, `homopolymer`, `direpeat`, `trirepeat`,
  and `quadrepeat`, are computed using the VCF fields `RU` and `RPA`.
  - If `RU` and `RPA` are not available, the features are computed manually
    using reference sequence 100 bp after the variant position.

| Name              | Source                         | Description                                                                                                                                         |
|:------------------|:-------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------|
| alt_ad            | field `AD`                     | Number of reads supporting the ALT-allele in the VCF file. This is the sum of both non-discordant and discordant supporting reads.                  |
| alt_ad2           | field `AD`                     | Number of reads supporting the other ALT-allele in the VCF file. This is the sum of both non-discordant and discordant supporting reads.            |
| alt_ad2_af        | fields `AD`,`DP`               | Allele frequency based on number of reads supporting the other ALT-allele in the VCF file; `alt_ad2` / `normal_dp`                                  |
| alt_ad_af         | fields `AD`,`DP`               | Allele frequency based on number of reads supporting the ALT-allele in the VCF file; `alt_ad` / `normal_dp`                                         |
| at_interest       | `--interest-regions` BED file  | Variant is in an interest region; `0`:False, `1`:True                                                                                               |
| direpeat          | reference or fields `RU`,`RPA` | Occurrence of 2-mer repeat after the variant's start position                                                                                       |
| gc_content        | reference                      | GC content of the reference sequence in a 200-bp window centered around the variant start position                                                  |
| genotype          | field `GT`                     | Genotype of the variant                                                                                                                             |
| hapcomp           | field `HAPCOMP`                | Edit distances of each ALT-allele's most common supporting haplotype from closest germline haplotype, excluding differences at the site in question |
| hapdom            | field `HAPDOM`                 | For each ALT-allele, fraction of read support that best fits the most-supported haplotype containing the allele                                     |
| mq_rank_sum       | INFO field `MQRankSum`         | GATK mapping-quality rank-sum annotation                                                                                                             |
| homopolymer       | reference or fields `RU`,`RPA` | Homopolymer length after the variant's start position                                                                                               |
| in_catalog        | `--catalog-vcf` VCF file       | Variant matches a record in the catalog VCF (chrom + position + trimmed REF/ALT); `0`:False, `1`:True                                                |
| normal_af         | field `AF`                     | Normal sample allele frequency from the VCF file                                                                                                    |
| normal_alt_ad     | field `AD`                     | Number of reads supporting the ALT-allele in normal sample only                                                                                     |
| normal_dp         | field `DP`                     | Depth of reads in normal sample only                                                                                                                |
| popaf             | `--pop-af-vcf` VCF file        | Population allele frequency for the ALT-allele                                                                                                      |
| post_2bp_context  | reference                      | 2-bp context after the variant's start position                                                                                                     |
| post_30bp_context | reference                      | 30-bp context after the variant's start position                                                                                                    |
| pre_2bp_context   | reference                      | 2-bp context before the variant's start position                                                                                                    |
| quadrepeat        | reference or fields `RU`,`RPA` | Occurrence of 4-mer repeat after the variant's start position                                                                                       |
| ref_ad            | field `AD`                     | Number of reads supporting the REF-allele in the VCF file. This is the sum of both non-discordant and discordant supporting reads.                  |
| ref_ad_af         | fields `AD`,`DP`               | Allele frequency based on number of reads supporting the REF-allele in the VCF file; `ref_ad` / `normal_dp`                                         |
| rpa_alt           | field `RPA`                    | Number of times tandem repeat unit is repeated for the ALT-allele                                                                                   |
| rpa_ref           | field `RPA`                    | Number of times tandem repeat unit is repeated for the REF-allele                                                                                   |
| sor               | INFO field `SOR`               | Symmetric Odds Ratio of 2x2 contingency table to detect strand bias                                                                                   |
| ru                | field `RU`                     | Tandem repeat unit (bases)                                                                                                                          |
| str               | field `STR`                    | Variant is a short tandem repeat; `0`:False, `1`:True                                                                                               |
| trirepeat         | reference or fields `RU`,`RPA` | Occurrence of 3-mer repeat after the variant's start position                                                                                       |
| tumor_af          | field `AF`                     | Tumor sample allele frequency from the VCF file                                                                                                     |
| tumor_alt_ad      | field `AD`                     | Number of reads supporting the ALT-allele in tumor sample only                                                                                      |
| tumor_dp          | field `DP`                     | Depth of reads in tumor sample only                                                                                                                 |
| uniq_3mers        | reference                      | Number of unique 3-mers in `post_30bp_context`                                                                                                      |
| uniq_4mers        | reference                      | Number of unique 4-mers in `post_30bp_context`                                                                                                      |
| uniq_5mers        | reference                      | Number of unique 5-mers in `post_30bp_context`                                                                                                      |
| uniq_6mers        | reference                      | Number of unique 6-mers in `post_30bp_context`                                                                                                      |
| variant_density   |                                | Number of variants within a 201-bp window (100 bp up/downstream) at the variant site                                                                |
| vcf_tn_af_ratio   | field `AF`                     | Ratio of tumor sample AF to sum of tumor and normal sample AF; `tumor_af` / (`tumor_af` + `normal_af`)                                              |
| vcf_variant_gq    | field `GQ`                     | Variant genotype quality                                                                                                                            |
| vcf_variant_qual  | `QUAL` column                  | Variant quality                                                                                                                                     |

Back to [compute_bam_features output](usage-compute-bam-features.md#output)
and [compute_vcf_features output](usage-compute-vcf-features.md#output).
