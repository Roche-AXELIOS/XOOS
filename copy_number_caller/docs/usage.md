# Usage

Use this page to choose the right Copy Number Caller command.
Follow the links below to navigate to the full command guides.

## Choosing a command

Use `germline-wgs` to call copy number variants from a germline normal (non-tumor) WGS sample.
Use `tumor-normal-wgs` to characterize allele-specific somatic CNVs from a matched tumor-normal pair.
Use `tumor-only-wgs` to characterize allele-specific somatic CNVs from a tumor sample when no matched normal is available.
Use `predict-somatic-cnv` to re-call copy numbers from existing `tumor-normal-wgs` or `tumor-only-wgs` outputs under a different purity/ploidy solution, without re-running the full pipeline.

Options shared across all commands are described in the [Appendix](appendix.md):
[Genome mappability](appendix.md#genome-mappability),
[Seed segments](appendix.md#seed-segments),
and [Allosome handling](appendix.md#allosome-handling).
See [Simplex and duplex base handling](../README.md#simplex-and-duplex-base-handling) in the overview for how base type is treated.

## Call germline copy numbers

Call copy number variants from a germline normal WGS sample.
Input: **BAM + reference** -> Output: **germline SEG, VCF, and visualization outputs**.

```bash
copy_number_caller germline-wgs \
        --bam INPUT.bam \
        --reference REFERENCE.fa \
        --output-dir OUTPUT_DIR
```

To run germline copy number calling, go to [Call germline copy numbers](usage-germline-wgs.md).

## Call tumor-normal somatic copy numbers

Characterize allele-specific somatic CNVs from a matched tumor-normal pair.
Input: **tumor BAM + normal BAM + VCF + reference** -> Output: **somatic SEG, VCF, purity/ploidy grid, and visualization outputs**.

```bash
copy_number_caller tumor-normal-wgs \
        --tumor-bam TUMOR.bam \
        --normal-bam NORMAL.bam \
        --vcf VARIANTS.vcf.gz \
        --reference REFERENCE.fa \
        --output-dir OUTPUT_DIR
```

To run somatic tumor-normal calling, go to [Call tumor-normal somatic copy numbers](usage-tumor-normal-wgs.md).

## Call tumor-only somatic copy numbers

Characterize allele-specific somatic CNVs from a tumor sample without a matched normal.
Input: **tumor BAM + tumor VCF + population VCF + reference** -> Output: **somatic SEG, VCF, purity/ploidy grid, and visualization outputs**.

```bash
copy_number_caller tumor-only-wgs \
        --tumor-bam TUMOR.bam \
        --vcf VARIANTS.vcf.gz \
        --population-vcf POPULATION.vcf.gz \
        --reference REFERENCE.fa \
        --output-dir OUTPUT_DIR
```

To run somatic tumor-only calling, go to [Call tumor-only somatic copy numbers](usage-tumor-only-wgs.md).

## Re-call somatic copy numbers

Re-call copy numbers from existing `tumor-normal-wgs` or `tumor-only-wgs` outputs under a user-specified purity/ploidy solution.
Input: **segments + log-ratios + BAFs from either workflow + reference + tumor purity + tumor ploidy + sex** -> Output: **somatic SEG and VCF**.

For inputs from `tumor-normal-wgs`, run:

```bash
copy_number_caller predict-somatic-cnv \
        --segments SEGMENTS.seg \
        --log-ratios LOG_RATIOS.bed \
        --bafs BAFS.bed \
        --reference REFERENCE.fa \
        --tumor-purity PURITY \
        --tumor-ploidy PLOIDY \
        --sex SEX \
        --output-dir OUTPUT_DIR
```

For inputs from `tumor-only-wgs`, run the command with `--tumor-only`. The mode is not inferred from the input files:

```bash
copy_number_caller predict-somatic-cnv \
        --segments SEGMENTS.seg \
        --log-ratios LOG_RATIOS.bed \
        --bafs BAFS.bed \
        --reference REFERENCE.fa \
        --tumor-purity PURITY \
        --tumor-ploidy PLOIDY \
        --sex SEX \
        --tumor-only \
        --output-dir OUTPUT_DIR
```

To re-call somatic CNVs, go to [Re-call somatic copy numbers](usage-predict-somatic-cnv.md).
