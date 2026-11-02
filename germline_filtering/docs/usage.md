# Usage

Germline Variant Filtering provides two commands: `train` builds classification models from labeled data, and
`predict` applies them to one sample. Choose the command based on whether you need to build models or classify
variants with existing models.

SNVs and indels are handled by separate models: `train` fits one model per variant type, and `predict` routes each variant to the model matching its type.

## Choosing a command

Use `train` when you have a labeled, multi-sample cohort and need to train models for a specific assay or capture method.
Custom training on WES or target-enrichment data is technically supported, but the module has not been tested on those sample types; training does not establish suitability for their use.
Use `predict` when you need to classify variants in one tumor-only sample with the shipped models or models
created by `train`.

## Train a model

Builds LightGBM classifiers from labeled multi-sample data described by a JSON config.
Input: **[labeled VCFs and SEG files via JSON](file-formats.md)** →
Output: **[SNV and indel model files](file-formats.md#model-files)**.

```shell
germline_filtering train --input training_samples.json --output-prefix my_experiment
```

To build models, go to [Train a model](usage-train.md).

## Apply a model to a sample

Applies trained models to one sample and produces a filtered, germline/somatic-tagged variant set for downstream somatic analysis.
Input: **[tumor-only VCF](file-formats.md#vcf-file) and
[copy-number SEG](file-formats.md#seg-file)** →
Output: **[tagged VCF](file-formats.md#tagged-vcf) and
[prediction CSV](file-formats.md#prediction-csv)**.

```shell
germline_filtering predict --vcf sample.vcf --seg sample.seg --output-prefix analysis_run
```

To classify variants in a sample, go to [Apply a model to a sample](usage-predict.md).

For input and output schemas, see [File formats](file-formats.md). For run summaries, see
[Metrics](metrics.md). If a run fails or produces unexpected counts, see [Troubleshooting](troubleshooting.md).
