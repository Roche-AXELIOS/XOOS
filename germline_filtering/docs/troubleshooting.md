# Troubleshooting

This page lists common issues and how to resolve them.
For the file requirements referenced below, see the [File formats](file-formats.md) page; for the CLI flags, see the [Train a model](usage-train.md) and [Apply a model to a sample](usage-predict.md) pages.

## Missing `##tumor_sample` header

**Symptom:** For a multi-sample VCF, the tool exits during VCF parsing with `Cannot determine the tumor sample column`. For a single-sample VCF, the tool logs a warning that the `##tumor_sample` value does not match the only sample column and continues.

**Cause:** The VCF has no `##tumor_sample=<sample_id>` header line, or the ID does not match a sample column in the `#CHROM` line (for example, after the sample column was renamed during reheadering).
A single-sample VCF falls back to its only sample column; a multi-sample VCF cannot be resolved unambiguously.

**Resolution:** Add or correct the header line so that `<sample_id>` exactly matches the tumor sample column name.
This is required for multi-sample VCFs and removes the warning for single-sample VCFs.
See [Required header lines](file-formats.md#required-header-lines).

## Depth filter drops too many variants

**Symptom:** `vcf_filtered_by_depth` in the QC metrics is high and `features_kept` (or `predictions_total`) is far lower than the number of input records.

**Cause:** Records with tumor `DP` below `--min-tumor-depth` (default 20) are dropped before scoring.

**Resolution:** Lower `--min-tumor-depth` if it is too aggressive for your data. Confirm the tumor sample declared in `##tumor_sample` actually carries `FORMAT/DP` values.

## Variants skipped for unusable FILTER values

**Symptom:** In `train`, `vcf_filtered_by_other` is high and few records reach training.

**Cause:** Only `PASS`, `germline`, and `germline;normal_artifact` are usable labels; everything else is skipped.

**Resolution:** Assign the correct germline/somatic `FILTER` labels to training records. For prediction, callable variants should carry `FILTER=PASS`. See [FILTER conventions (labels)](file-formats.md#filter-conventions-labels).

## Variants have missing POPAF

**Symptom:** Many inputs lack an `INFO/POPAF` value and predictions rely more heavily on other features.

**Cause:** Missing POPAF is NaN-encoded and routed by the model; it no longer prevents feature generation.

**Resolution:** Population-frequency annotation is recommended because it is highly discriminative, but it is not required. To add it, use a command such as `bcftools annotate -a gnomad.vcf.gz -c INFO/POPAF`. See [How to create these VCF files](file-formats.md#how-to-create-these-vcf-files).

## NaN feature values drop variants

**Symptom:** `features_skipped_na` is non-zero and fewer rows reach scoring than expected.

**Cause:** One or more feature cells could not be computed, for example a missing SEG value for the segment covering the variant.

**Resolution:** Confirm the SEG file is genome-wide and well formed, with valid `Purity`, `Ploidy`, and copy-number columns. Variants outside all SEG segments use a diploid fallback during feature generation, so persistent NA usually points to malformed SEG rows. In `predict`, SEG-gap records are emitted as `nocall_segment_missing` and excluded from both final call sets. See [SEG file](file-formats.md#seg-file).

## Reference-genome or coordinate mismatch

**Symptom:** Copy-number features look wrong, or `features_cnv_gap_fallback` is unexpectedly high.

**Cause:** The VCF and SEG files use different reference builds or coordinate systems.

**Resolution:** Ensure both files use the same reference genome (GRCh38 is tested) and the same coordinate convention.

## Model fails to load or produces nonsensical scores

**Symptom:** `predict` errors while loading a model, or every variant scores near an extreme.

**Cause:** The model was trained on an incompatible feature layout.

**Resolution:** Retrain with the current version, or use a model produced by the current version. See [Model files](file-formats.md#model-files).

## Overrides are not tagging any variants

**Symptom:** `predictions_high_af` or `predictions_common_germline` is always 0.

**Cause:** The overrides are disabled by default; `--max-somatic-af` and `--common-germline-popaf` both default to `1.0`, which turns the respective override off.

**Resolution:** Set `--max-somatic-af` below `1.0` to enable the AF override, and `--common-germline-popaf` (and, for indels, `--indel-common-germline-popaf`) below `1.0` to enable the POPAF override. See [Classification overrides](appendix.md#classification-overrides).

## Prediction is slow on large VCFs

**Symptom:** `predict` takes a long time on whole-genome inputs.

**Cause:** The default `--threads 1` uses the single-threaded path.

**Resolution:** Increase `--threads` (use `--threads 0` for all hardware threads). The VCF is partitioned into genomic regions scored in parallel and concatenated in order; decompressed VCF content is identical regardless of thread count.

## Training log shows overfitting

**Symptom:** Validation `binary_logloss` rises while training loss keeps falling.

**Cause:** The model is fitting noise, typically with too little labeled data per variant type.

**Resolution:** Early stopping is on by default (`--early-stopping-rounds 50`) and truncates the model to the best iteration. If needed, lower `--num-iterations` or increase `--early-stopping-rounds`, and confirm you have enough labeled data for each of the separately-trained SNV and indel models.
