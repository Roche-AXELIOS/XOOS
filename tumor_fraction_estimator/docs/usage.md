# Usage

Use this page to choose the right Tumor Fraction Estimator subcommand.
Follow the links below to navigate to the full command guides.

## Choosing between `tumor-fraction` and `contamination`

Both subcommands share the same analytical engine.
They count reference and alternate alleles at a set of query variant sites, apply the same families of filters, and estimate a fraction as `2 * mean_vaf`.
They differ in the input variant set, the default filter thresholds, and the output labels.

Use `tumor-fraction` to estimate the fraction of circulating tumor DNA from a somatic variant VCF that has a matched tumor and normal sample.
Use `contamination` to estimate cross-sample contamination from population SNP sites.
The bundled 1000 Genomes population SNP VCF is used by default, so `contamination` can run without a supplied VCF.

## Estimate tumor fraction

Estimate tumor fraction from somatic variant sites in cfDNA sequencing data.
Input: **BAM + somatic VCF + reference** -> Output: **`tumor_fraction_estimation_results.tsv` and `per_site_counts.tsv`**.

```shell
tumor_fraction_estimator tumor-fraction \
  --bam INPUT.bam \
  --query-vcf SOMATIC.vcf \
  --tumor-sample-name TUMOR \
  --reference REFERENCE.fa \
  --included-regions-bed REGIONS.bed \
  --output-dir OUTPUT_DIR
```

To run tumor fraction estimation, go to [Estimate tumor fraction](usage-tumor-fraction.md).

## Estimate contamination

Estimate cross-sample contamination from population SNP sites.
Input: **BAM + reference (+ optional VCF)** -> Output: **`contamination_estimation_results.tsv` and `per_site_counts.tsv`**.

```shell
tumor_fraction_estimator contamination \
  --bam INPUT.bam \
  --reference REFERENCE.fa \
  --included-regions-bed REGIONS.bed \
  --output-dir OUTPUT_DIR
```

To run contamination estimation, go to [Estimate contamination](usage-contamination.md).
