# Metrics

Both subcommands write a per-sample QC metrics file (enabled by default; pass `--no-qc-metrics` to skip).
The metrics summarize the input VCF and SEG, the feature-processing pipeline, and — for `predict` — the classification and tagged-VCF outputs.

Prediction performance metrics are opt-in. Pass `--performance-metrics` to write
`<output-prefix-parent>/metric/<sample>.performance.tsv`; the option is disabled by default and, since it is
infrequently needed, only shown under `predict --help-advanced`. This requires a **labeled validation
VCF** with both `FILTER=PASS` and `FILTER=germline[;normal_artifact]` records -- a raw Tumor-Only VCF
(all `PASS`) fails this requirement; see [Prediction performance TSV](#prediction-performance-tsv). This
failure is logged as an error and skips only the performance metrics file: the rest of predict's
outputs, including the QC metrics TSV above, are still written and the run still exits successfully.

## QC metrics TSV

Applies to: `train`, `predict`.
Default path: `<output-prefix-parent>/metric/<sample>.tsv`.

Per-sample QC metrics are written to a `metric/` directory next to the output prefix, one file per sample named
`<sample>.tsv`. For `train`, the sample names come from the JSON configuration. For `predict`, the sample name is
the final path component of `--output-prefix`. Each file is a two-column `metric<TAB>value` TSV.

### Input VCF counts

| Metric | Definition |
|:-------|:-----------|
| vcf_total_records | Total records read from the input VCF. |
| vcf_pass_variants | Records with `FILTER=PASS`. In `train`, these are labeled somatic variants. In `predict`, they are the callable variants submitted for germline/somatic classification; their input FILTER value does not indicate the final classification. |
| vcf_germline_variants | Records with a germline `FILTER` label. |
| vcf_filtered_by_depth | Records dropped for tumor `DP` below `--min-tumor-depth`. |
| vcf_filtered_by_other | Records skipped because the `FILTER` value is not a usable label. |
| vcf_phased_variants_skipped_invalid_mnvtag | Phased compound records skipped because `INFO/MNVTAG` is missing or its first variant cannot be parsed. |
| vcf_snv_variants | SNVs among kept records. |
| vcf_indel_variants | indels among kept records. |
| vcf_sex_chromosome_variants | Kept records on `chrX`/`chrY`. |
| vcf_variants_\<chromosome\> | Kept-record count per chromosome. |

### Input SEG summary

| Metric | Definition |
|:-------|:-----------|
| seg_num_segments | Number of segments in the SEG file. |
| seg_segments_in_allosome | Segments on sex chromosomes. |
| seg_mean_total_copy_number | Mean total copy number across segments. |
| seg_mean_minor_copy_number | Mean minor allele copy number across segments. |
| seg_loh_fraction | Fraction of the assessable autosomal genome (by length) with loss of heterozygosity. Allosomes and segments without allelic evidence are excluded from both sides of the ratio. |
| seg_segments_without_baf | Segments carrying no allelic (BAF) evidence. A high count means the LOH features are unknown for much of the genome and the copy-number caller had thin SNP support. |
| seg_purity | Tumor purity from the SEG file. |
| seg_ploidy | Tumor ploidy from the SEG file. |

### Feature-processing counts

| Metric | Definition |
|:-------|:-----------|
| features_before_deduplication | Feature rows before deduplication. |
| features_duplicates_removed | Duplicate feature rows removed within each sample. |
| features_after_deduplication | Feature rows after deduplication. |
| features_cnv_gap_fallback | Variants using the diploid CNV-gap fallback. |
| features_total_rows | Feature rows entering the pre-scoring filter. |
| features_skipped_sex_chromosome | Rows dropped as sex-chromosome variants. |
| features_indel_scored | indel rows flowing to scoring. |
| features_skipped_na | Rows dropped for NA feature values. |
| features_kept | Rows kept for scoring. |

### Prediction counts (`predict` only)

| Metric | Definition |
|:-------|:-----------|
| predictions_total | Total scored variants. |
| predictions_germline | Variants classified `germline` by the model. |
| predictions_high_af | Variants reclassified `germline` by the AF override. |
| predictions_common_germline | Variants reclassified `germline` by the POPAF override. |
| predictions_somatic | Variants kept somatic. |

### Tagged-VCF counts (`predict` only)

These are the **final, no-call-aware** output counts. They can differ from the `predictions_*` counts above because a no-call/ambiguous status overrides the model germline/somatic decision (`predictions_*` are the model decisions; `tagged_vcf_*` are what the tagged VCF actually emits).

| Metric | Definition |
|:-------|:-----------|
| tagged_vcf_total_records | Total records written to the tagged VCF. |
| tagged_vcf_germline_tagged | Records tagged `germline` by the model. |
| tagged_vcf_high_af_tagged | Records tagged `high_af` by the AF override. |
| tagged_vcf_common_germline_tagged | Records tagged `common_germline` by the POPAF override. |
| tagged_vcf_somatic_tagged | Records kept somatic (`PASS`). |
| tagged_vcf_nocall_cyp2d6 | Records no-called for the CYP2D6/CYP2D7 locus. |
| tagged_vcf_nocall_hla | Records no-called for the MHC/HLA region. |
| tagged_vcf_nocall_purity | Records no-called because sample purity exceeds 0.95. |
| tagged_vcf_nocall_segment_missing | Scored records no-called because no copy-number segment covers them. |
| tagged_vcf_ambiguous_cna | Scored records marked ambiguous because the covered segment's copy number is untrustworthy. |
| tagged_vcf_nocall_coverage | Unscored records no-called because tumor `DP` is below `--min-tumor-depth`. |
| tagged_vcf_unknown | Unscored records not caught above (sex chromosome / NaN features); `FILTER` unchanged. |
| tagged_vcf_dropped | Records dropped under `--no-emit-all`. |

## Prediction performance TSV

**Class convention**: germline is the **positive** class and somatic (`FILTER=PASS`) is the
**negative** class, matching the truth-label convention used elsewhere (a germline `FILTER` label is
the condition of interest; `PASS` is the default/negative call).

**Input requirement**: the input VCF must be a **labeled validation VCF** containing both
`FILTER=PASS` (truth-negative) and `FILTER=germline[;normal_artifact]` (truth-positive) records.
A raw Tumor-Only (TO) VCF -- where every real variant carries `FILTER=PASS` because TO input does not
distinguish germline from somatic, which is the whole reason `germline_filtering` exists -- cannot be
used: predict only scores `PASS` candidates, so a real germline variant it correctly flags would be
counted as a false positive against a truth label that is wrong for that variant, producing meaningless
rates (for example, a false-positive rate near 1.0 that reflects nothing but the TO VCF's mixed truth
population). The evaluation requires at least one record of each label; if that requirement isn't met,
predict logs an error and skips writing this file, but the rest of the run (including the QC metrics
TSV) still completes and the process still exits successfully.

Predict's feature extraction is `PASS`-only: a germline-labeled (truth-positive) input record is
never scored by the model. This file therefore can only measure specificity/false-positive behavior
on truth-negative (`FILTER=PASS`) records predict did score -- how many were correctly kept somatic
(true negative) vs. incorrectly reclassified germline (false positive); an AF/POPAF override status
counts as a false positive, since it is a final germline reclassification. It does **not** measure
germline-detection sensitivity/recall (true positives, false negatives): no germline-labeled
(truth-positive) record ever reaches the model, so those records cannot appear in the overall/SNV/indel
counts below and are reported only in the top-level `evaluation_germline_truth_records` summary count,
for transparency. No-call and `unknown` statuses are excluded from the somatic-retention counts and
reported in the corresponding `evaluation_skipped_*` counters. SNV and indel counts are reported
separately as well as in the `overall_*` metrics.

| Metric | Definition |
|:-------|:-----------|
| evaluation_total_input_records | Total records read from the input VCF. |
| evaluation_labeled_records | Records with a usable truth label (`FILTER=PASS` or `FILTER=germline[;normal_artifact]`). |
| evaluation_somatic_truth_records | Truth-negative (`FILTER=PASS`) records; the only class predict scores and evaluates below. |
| evaluation_germline_truth_records | Truth-positive (germline) records. Informational only: never scored by predict, so never part of the counts below. |
| evaluation_records | Truth-negative records that were successfully evaluated (i.e. counted in `overall_somatic_truth_evaluated`). |
| evaluation_skipped_unlabeled | Records with no usable `FILTER` label. |
| evaluation_skipped_not_emitted | Truth-negative records absent from the tagged VCF (dropped by `--no-emit-all`, though `PASS` records are never dropped by it). |
| evaluation_skipped_no_call | Truth-negative records with a no-call status. |
| evaluation_skipped_unknown | Truth-negative records with an `unknown` status (unscored: sex chromosome / NaN features). |

For each of `overall`, `snv`, and `indel`:

| Metric | Definition |
|:-------|:-----------|
| \<prefix\>_somatic_truth_evaluated | Truth-negative records evaluated (true negatives + false positives). |
| \<prefix\>_correctly_retained_somatic | True negatives: truth-negative records the final classification kept somatic. |
| \<prefix\>_false_germline_calls | False positives: truth-negative records the final classification reclassified as germline. |
| \<prefix\>_false_germline_rate | False positive rate: `false_germline_calls / somatic_truth_evaluated`. `NA` when the denominator is zero. |
| \<prefix\>_somatic_retention_rate | Specificity (true negative rate): `correctly_retained_somatic / somatic_truth_evaluated`. `NA` when the denominator is zero. |
