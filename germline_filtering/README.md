# Germline Variant Filtering

> **Warning:** Germline Variant Filtering is for research use only.

## Guide contents

- [Overview](#overview): what the tool does, its model-based classification, example applications, and system requirements.
- [Usage](docs/usage.md): running each subcommand.
  - [Train a model](docs/usage-train.md): building SNV/indel models from labeled data, with standard `train` CLI options.
  - [Apply a model to a sample](docs/usage-predict.md): running `predict` on a sample, with standard `predict` CLI options.
- [File formats](docs/file-formats.md): input JSON/VCF/SEG and output model/CSV/tagged-VCF layouts.
- [Metrics](docs/metrics.md): the per-sample QC metrics TSV and its fields.
- [Troubleshooting](docs/troubleshooting.md): common issues and resolutions.
- [Appendix](docs/appendix.md): concepts and terminology, feature engineering, classification overrides, and training internals.
- [Release notes](CHANGELOG.md): changes by release.

## Overview

Germline Variant Filtering classifies variants in tumor-only samples as likely germline or somatic for research
analysis. Removing likely germline variants helps downstream analyses focus on the somatic variant set.

This tool uses machine learning (LightGBM gradient boosting) to classify variants as germline or somatic.
Each variant becomes a 31-feature vector combining VCF evidence (i.e. allele fraction, read depth, variant clustering, population frequency) with copy-number context from the SEG file (i.e. copy number, purity, ploidy, loss of heterozygosity and allelic imbalance at the locus, the sample's genome-wide LOH burden, and a purity-normalized allele fraction).
SNVs and indels are scored by separate models, and two deterministic overrides can reclassify a model-somatic call as germline outside the model.
It exposes two subcommands:

- **`train`** — builds classification models from labeled training data described by a JSON config.
- **`predict`** — applies trained models to a single sample's VCF/SEG inputs to classify variants and produce a filtered, tagged variant set for downstream somatic analysis.

At a glance: **tumor-only VCF + copy-number SEG** →
**[germline/somatic-tagged VCF](docs/file-formats.md#tagged-vcf)** (`predict`).

The tool has two workflows. Use [Train a model](docs/usage-train.md) to build SNV and indel models from labeled
cohort data. Use [Apply a model to a sample](docs/usage-predict.md) to apply trained or shipped models to one sample.
See the [Usage](docs/usage.md) page for command-selection guidance.

### Limitation

Germline Variant Filtering is intended for tumor-only variant-calling workflows and is for research use only.
It has been tested only on whole-genome sequencing (WGS) FFPE samples.
Use with whole-exome sequencing (WES) or target-enrichment samples is outside the tested scope and is not
recommended.
Training a custom model with WES or target-enrichment data is technically supported, but does not establish that
the module or resulting model is suitable for those sample types; independent evaluation is required.

Prediction requires a genome-wide copy-number segmentation (SEG) file for the same sample, produced by the Copy
Number Caller.
Variants without a covering segment, or with insufficient or ambiguous copy-number evidence, may be reported as
no-calls and excluded from the germline and somatic sets.
The tool classifies autosomal variants only; sex-chromosome variants are excluded.

### Model-based classification

`train` builds separate LightGBM models for SNVs and indels from labeled VCF and copy-number SEG files.
For each variant, the model combines variant-level evidence—allele fraction, read depth, population frequency,
and nearby common-variant density—with sample and locus-specific copy-number context. That context includes
purity, ploidy, total and allele-specific copy number, loss of heterozygosity, allelic imbalance, and expected
allele fractions under germline and somatic hypotheses.

The tool calculates 31 features per variant and uses 23 of them for model training and prediction. The remaining
features are retained for output and bookkeeping. See [Feature engineering](docs/appendix.md#feature-engineering)
for the complete feature definitions and [Model features](docs/appendix.md#model-features) for the active model
feature set.

### Key features

- Machine learning classification using LightGBM with 23 model features (from a 31-feature physical matrix).
- Separate SNV and indel models; each variant is routed to the matching model.
- Integration of tumor purity, ploidy, and copy-number alterations from the SEG file.
- Loss-of-heterozygosity context, both per locus and as a genome-wide burden, so germline variants whose allele fraction is skewed by allelic loss are not mistaken for somatic calls. See [LOH features](docs/appendix.md#loh-features).
- Diploid copy-number fallback (CN=2, minor allele=1) for feature generation when variants fall outside SEG
  segments. In `predict`, these records are emitted as `nocall_segment_missing` rather than being included in the
  germline or somatic sets.
- Phased compound records (`PASS`, `len(REF) > 1` and `len(ALT) > 1`) are routed to the SNV or indel model
  based on the first variant listed in `INFO/MNVTAG`, not the compound record's own REF/ALT length. Records
  without a usable first `MNVTAG` variant are skipped.
- Optional deterministic overrides, disabled by default, are available for specialized analyses. See
  [Classification overrides](docs/appendix.md#classification-overrides).
- No-call statuses that abstain on untrustworthy loci (CYP2D6/HLA regions, near-pure samples, missing/ambiguous copy number, and below-depth records), excluding them from both the germline and somatic sets. The final status is written to `INFO/GERMLINE_STATUS`. See [No-call statuses](docs/appendix.md#no-call-statuses).
- Multithreaded prediction with output identical regardless of thread count.
- Per-sample QC metrics (input, feature-processing, prediction, and tagged-VCF counts).
- Sex-chromosome filtering (autosomes only).

### Example applications

Each example includes its data type, objective, command, and detailed usage guide.

#### Filter germline variants from a tumor-only sample

- **Data type:** Tumor-only VCF with a matched copy-number SEG.
- **Objective:** Exclude likely germline variants from a downstream research analysis of somatic calls, using the models shipped with the tool.
- **Usage:** To run this application, go to [Apply a model to a sample](docs/usage-predict.md).
- **Command:** `germline_filtering predict --vcf sample.vcf --seg sample.seg --output-prefix analysis_run`

#### Train a capture-specific classifier

- **Data type:** Labeled multi-sample cohort (VCF + SEG per sample) described by a JSON config.
- **Objective:** Train SNV and indel models on a cohort from a specific capture method. Training on WES or target-enrichment data is technically supported, but is outside the module's tested scope and does not establish suitability for those sample types; independent evaluation is required.
- **Usage:** To run this application, go to [Train a model](docs/usage-train.md).
- **Command:** `germline_filtering train --input training_samples.json --output-prefix my_experiment`

### System requirements

|          | Requirement |
|:---------|:------------|
| CPU      | Modern processor with multiple cores (prediction is multithreaded). |
| Memory   | At least 16 GiB. |
| Software | Docker or an equivalent containerization environment. |
