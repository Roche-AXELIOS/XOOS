# Generating candidate calls

The upstream caller must emit both the candidate VCF and the realigned BAM that
SVC uses for feature computation.
The following examples provide a GATK command and an NVIDIA Parabricks command
for each workflow.
The NVIDIA Parabricks examples use syntax supported by Parabricks 4.7.1.

## Germline whole-genome

The supported library types are SBX-Simplex (SBX-S), SBX-Duplex (SBX-D), and
SBX-Fast.

### Duplex data (SBX-Duplex (SBX-D) and SBX-Fast adapters)

Use GATK HaplotypeCaller or NVIDIA Parabricks `haplotypecaller` to generate
inputs for duplex `germline-wgs` runs.
See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs)
for the corresponding `--library-type` values.

The input *ALIGNMENT.bam* must be a duplicate-marked BAM, for example, the
output from Read Collapser `markdup`, and not an inter-consensus BAM.

Run GATK as follows:

```bash
gatk HaplotypeCaller \
    -I ALIGNMENT.bam \
    -R REFERENCE.fa \
    --intervals TARGET.bed \
    -O GATK_OUTPUT.vcf.gz \
    -bamout GATK_OUTPUT.bam \
    --bam-writer-type CALLED_HAPLOTYPES \
    -A AssemblyComplexity \
    -A TandemRepeat \
    -RF MappingQualityReadFilter \
    --activeregion-alt-multiplier 5 \
    --adaptive-pruning true \
    --enable-dynamic-read-disqualification-for-genotyping true \
    --mapping-quality-threshold-for-genotyping 1 \
    --minimum-mapping-quality 1 \
    --min-base-quality-score 6 \
    --native-pair-hmm-threads 4 \
    --smith-waterman FASTEST_AVAILABLE \
    --max-mnp-distance 0 \
    --create-output-bam-index true \
    --create-output-variant-index true \
    --tmp-dir /PATH/TO/TEMP/DIRECTORY
```

Key parameters for germline variant filtering:

- `-A AssemblyComplexity` and `-A TandemRepeat` produce additional VCF
  annotations for machine learning features in SVC.
- `--min-base-quality-score 6` skips discordant bases in duplex consensus
  reads.
- `--activeregion-alt-multiplier 5` increases sensitivity.
  Omit this option for SBX-Simplex (SBX-S) data.
- `--max-mnp-distance 0` splits multi-nucleotide variants (MNVs) into SNVs.
  - This is required to properly filter and re-genotype SNVs.

See the [GATK HaplotypeCaller user guide](https://gatk.broadinstitute.org/hc/en-us/articles/9570334998171-HaplotypeCaller)
for other parameters.

Run NVIDIA Parabricks as follows:

```bash
pbrun haplotypecaller \
    --ref REFERENCE.fa \
    --in-bam ALIGNMENT.bam \
    --interval-file TARGET.bed \
    --out-variants PARABRICKS_OUTPUT.vcf.gz \
    --htvc-bam-output PARABRICKS_OUTPUT.bam \
    --min-base-quality-score 6 \
    --num-htvc-threads THREADS \
    --mapping-quality-threshold-for-genotyping 1 \
    --minimum-mapping-quality 1 \
    --adaptive-pruning \
    --enable-dynamic-read-disqualification-for-genotyping \
    --filter-reads-too-long \
    --haplotypecaller-options '-A AssemblyComplexity -A TandemRepeat --activeregion-alt-multiplier 5'

samtools index PARABRICKS_OUTPUT.bam
bcftools index --tbi PARABRICKS_OUTPUT.vcf.gz
```

The current NVIDIA Parabricks HaplotypeCaller documentation does not list an
option equivalent to the explicit GATK `--max-mnp-distance 0` setting.
Validate the caller outputs before treating them as equivalent.

See the [NVIDIA Parabricks HaplotypeCaller user guide](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/haplotypecaller)
for other parameters.

### Simplex data (SBX-Simplex (SBX-S), YS adapters)

Use GATK HaplotypeCaller or NVIDIA Parabricks `haplotypecaller` to generate
inputs for simplex `germline-wgs` runs.
See [Library types for germline WGS](appendix.md#library-types-for-germline-wgs)
for the corresponding `--library-type` value.

The input *ALIGNMENT.bam* must be a duplicate-marked BAM, for example, the
output from Read Collapser `markdup`.

Run GATK as follows:

```bash
gatk HaplotypeCaller \
    -I ALIGNMENT.bam \
    -R REFERENCE.fa \
    --intervals TARGET.bed \
    -O GATK_OUTPUT.vcf.gz \
    -bamout GATK_OUTPUT.bam \
    -A AssemblyComplexity \
    -A TandemRepeat \
    -RF MappingQualityReadFilter \
    --adaptive-pruning true \
    --enable-dynamic-read-disqualification-for-genotyping true \
    --mapping-quality-threshold-for-genotyping 1 \
    --minimum-mapping-quality 1 \
    --min-base-quality-score 6 \
    --native-pair-hmm-threads 16 \
    --smith-waterman FASTEST_AVAILABLE \
    --create-output-bam-index true \
    --create-output-variant-index true
```

See the [GATK HaplotypeCaller user guide](https://gatk.broadinstitute.org/hc/en-us/articles/9570334998171-HaplotypeCaller)
for other parameters.

Run NVIDIA Parabricks as follows:

```bash
pbrun haplotypecaller \
    --ref REFERENCE.fa \
    --in-bam ALIGNMENT.bam \
    --out-variants PARABRICKS_OUTPUT.vcf.gz \
    --htvc-bam-output PARABRICKS_OUTPUT.bam \
    -G StandardAnnotation \
    -G StandardHCAnnotation \
    --num-gpus GPUS \
    --run-partition \
    --num-htvc-threads 6 \
    --min-base-quality-score 6 \
    --haplotypecaller-options '-A TandemRepeat -A AssemblyComplexity --adaptive-pruning --mapping-quality-threshold-for-genotyping 1 --minimum-mapping-quality 1 --enable-dynamic-read-disqualification-for-genotyping --filter-reads-too-long'

samtools index PARABRICKS_OUTPUT.bam
bcftools index --tbi PARABRICKS_OUTPUT.vcf.gz
```

NVIDIA Parabricks HaplotypeCaller cannot process reads longer than 500 bp.
Specify `--filter-reads-too-long` through `--haplotypecaller-options` to filter
these reads and prevent errors.

See the [NVIDIA Parabricks HaplotypeCaller user guide](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/haplotypecaller)
for other parameters.

## Germline mitochondria

Use GATK Mutect2 or NVIDIA Parabricks `mutectcaller` in mitochondria mode to
generate candidate mitochondrial variants.
The input BAM file to Mutect2/mutectcaller is either the output from [mitochondrial bam_post_processing steps](../../bam_post_process/README.md) run on an input Giraffe BAM or a BWA aligned BAM.
The germline mitochondria workflow currently supports duplex data only.
Adjust the `chrM` interval to match the mitochondrial contig name in the
reference.

Run GATK as follows:

```bash
gatk Mutect2 \
    -I INPUT.bam \
    -O GATK_OUTPUT.vcf.gz \
    -R REFERENCE.fa \
    --mitochondria-mode \
    -L chrM \
    --native-pair-hmm-threads 10 \
    --bamout GATK_OUTPUT.bam \
    --bam-writer-type CALLED_HAPLOTYPES \
    --max-mnp-distance 0 \
    --min-base-quality-score 23 \
    -RF MappingQualityReadFilter \
    --minimum-mapping-quality 1 \
    --create-output-bam-index true \
    --create-output-variant-index true
```

Key parameters for germline mitochondrial calling:

- `--mitochondria-mode`: enables the Mutect2 preset tuned for high-depth
  mitochondrial DNA.
- `-L chrM` restricts calling to the mitochondrial contig.
  Adjust the contig name to match the reference, for example, `chrM`.
- `--native-pair-hmm-threads 10` sets the number of PairHMM threads.
- `--bamout GATK_OUTPUT.bam`: emits the realigned BAM that SVC consumes for BAM
  feature extraction.
- `--max-mnp-distance 0`: disables merging of adjacent SNVs into multi-nucleotide
  variants, so each SNV is emitted as a separate record for SVC ingestion.
- `--min-base-quality-score 23`: minimum base quality for a base to be
  considered.
- `-RF MappingQualityReadFilter --minimum-mapping-quality 1`: drops reads with
  mapping quality `0`.

See the [GATK Mutect2 user guide](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
for other parameters.

Run NVIDIA Parabricks as follows:

```bash
pbrun mutectcaller \
    --ref REFERENCE.fa \
    --in-tumor-bam INPUT.bam \
    --tumor-name SAMPLE_NAME \
    --out-vcf PARABRICKS_OUTPUT.vcf.gz \
    --mutect-bam-output PARABRICKS_OUTPUT.bam \
    --mitochondria-mode \
    --interval chrM \
    --num-htvc-threads THREADS \
    --max-mnp-distance 0 \
    --min-base-quality-score 23 \
    --minimum-mapping-quality 1 \
    --filter-reads-too-long

samtools index PARABRICKS_OUTPUT.bam
bcftools index --tbi PARABRICKS_OUTPUT.vcf.gz
```

See the [NVIDIA Parabricks Mutectcaller user guide](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/mutectcaller)
for other parameters.

## Tumor-normal whole-genome

Use GATK Mutect2 4.3 or NVIDIA Parabricks `mutectcaller` 4.7.1 to generate
candidate variants from a matched tumor-normal pair.
Other GATK versions can produce missing or incompatible VCF annotations that
affect SVC filtering.

The input *TUMOR_ALIGNMENT.bam* and *NORMAL_ALIGNMENT.bam* files must be
duplicate-marked BAMs, for example, outputs from Read Collapser `markdup`, and
not inter-consensus BAMs.

Run GATK as follows:

```bash
gatk Mutect2 \
    -R REFERENCE.fa \
    -I TUMOR_ALIGNMENT.bam \
    -I NORMAL_ALIGNMENT.bam \
    --tumor TUMOR_SAMPLE_NAME \
    --normal NORMAL_SAMPLE_NAME \
    --output GATK_OUTPUT.vcf.gz \
    --bam-output GATK_OUTPUT.bam \
    --bam-writer-type CALLED_HAPLOTYPES \
    --native-pair-hmm-threads 4 \
    --min-base-quality-score 35 \
    --base-quality-score-threshold 35 \
    --minimum-mapping-quality 1 \
    --f1r2-median-mq 1 \
    --initial-tumor-lod -5 \
    --tumor-lod-to-emit -5 \
    --pruning-lod-threshold -9 \
    --max-reads-per-alignment-start 0 \
    --active-probability-threshold 0.00005 \
    --smith-waterman FASTEST_AVAILABLE \
    -A AssemblyComplexity \
    -A TandemRepeat \
    --normal-lod 2.2 \
    --min-pruning 2 \
    --max-mnp-distance 0 \
    --min-dangling-branch-length 0 \
    --allow-non-unique-kmers-in-ref true \
    --enable-dynamic-read-disqualification-for-genotyping true \
    --recover-all-dangling-branches true \
    --linked-de-bruijn-graph false \
    --pileup-detection false \
    --genotype-germline-sites true \
    --create-output-bam-index true \
    --create-output-variant-index true
```

Key parameters for somatic variant detection:

- `-A AssemblyComplexity` and `-A TandemRepeat` produce additional VCF
  annotations for machine learning features in SVC.
- `--min-base-quality-score 35` and `--base-quality-score-threshold 35` use only
  concordant bases in duplex consensus reads.
- `--max-mnp-distance 0` splits MNVs into SNVs.
  - This is required to properly filter and re-genotype SNVs.

See the [GATK Mutect2 user guide](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
for other parameters.

Run NVIDIA Parabricks as follows:

```bash
pbrun mutectcaller \
    --ref REFERENCE.fa \
    --in-tumor-bam TUMOR_ALIGNMENT.bam \
    --tumor-name TUMOR_SAMPLE_NAME \
    --in-normal-bam NORMAL_ALIGNMENT.bam \
    --normal-name NORMAL_SAMPLE_NAME \
    --mutect-bam-output PARABRICKS_OUTPUT.bam \
    --out-vcf PARABRICKS_OUTPUT.vcf.gz \
    --num-htvc-threads THREADS \
    --minimum-mapping-quality 1 \
    --min-base-quality-score 35 \
    --f1r2-median-mq 1 \
    --base-quality-score-threshold 35 \
    --allow-non-unique-kmers-in-ref \
    --enable-dynamic-read-disqualification-for-genotyping \
    --recover-all-dangling-branches \
    --mutectcaller-options '-A AssemblyComplexity -min-dangling-branch-length 0 -max-reads-per-alignment-start 0' \
    --filter-reads-too-long \
    --active-probability-threshold 0.00005 \
    --pruning-lod-threshold -9 \
    --tumor-lod-to-emit -5 \
    --initial-tumor-lod -5 \
    --max-mnp-distance 0 \
    --genotype-germline-sites \
    --normal-lod 2.2

samtools index PARABRICKS_OUTPUT.bam
bcftools index --tbi PARABRICKS_OUTPUT.vcf.gz
```

### Parabricks compatibility

SVC was tested with output from NVIDIA Parabricks `mutectcaller` 4.7.1.
The following differences apply to the Parabricks command:

- The current NVIDIA Parabricks Mutectcaller documentation lists
  `--pileup-detection`, but the command omits it to preserve the validated SVC
  recipe.
  Validate caller output before adding this option.
- The current documentation does not list `--min-pruning`, so the command uses
  the NVIDIA Parabricks default.
- `--recover-all-dangling-branches` is a flag without a value in Parabricks,
  whereas GATK takes a Boolean value.

See the [NVIDIA Parabricks Mutectcaller user guide](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/mutectcaller)
for other parameters.

## Tumor-only whole-genome

Use GATK Mutect2 4.3 or NVIDIA Parabricks `mutectcaller` 4.7.1 without a
matched normal to generate candidate variants for the `tumor-only-wgs`
workflow.
Other GATK versions can produce missing or incompatible VCF annotations that
affect SVC filtering.
The input BAM must be duplicate-marked and compatible with the reference used
for the Panel of Normals.

Run GATK as follows:

```bash
gatk Mutect2 \
    --input INPUT.bam \
    --reference REFERENCE.fa \
    --output GATK_OUTPUT.vcf.gz \
    --bamout GATK_OUTPUT.bam \
    --bam-writer-type CALLED_HAPLOTYPES \
    --native-pair-hmm-threads 4 \
    --min-base-quality-score 35 \
    --base-quality-score-threshold 35 \
    --minimum-mapping-quality 1 \
    --f1r2-median-mq 1 \
    --initial-tumor-lod -5 \
    --tumor-lod-to-emit -5 \
    --pruning-lod-threshold -9 \
    --max-reads-per-alignment-start 0 \
    --active-probability-threshold 0.00005 \
    --smith-waterman FASTEST_AVAILABLE \
    --max-mnp-distance 0 \
    -A AssemblyComplexity \
    -A TandemRepeat \
    --min-pruning 2 \
    --min-dangling-branch-length 0 \
    --allow-non-unique-kmers-in-ref true \
    --enable-dynamic-read-disqualification-for-genotyping true \
    --recover-all-dangling-branches true \
    --linked-de-bruijn-graph false \
    --pileup-detection false \
    --genotype-germline-sites true \
    --create-output-bam-index true \
    --create-output-variant-index true
```

Run NVIDIA Parabricks as follows:

```bash
pbrun mutectcaller \
    --ref REFERENCE.fa \
    --in-tumor-bam INPUT.bam \
    --tumor-name TUMOR_SAMPLE_NAME \
    --out-vcf PARABRICKS_OUTPUT.vcf.gz \
    --mutect-bam-output PARABRICKS_OUTPUT.bam \
    --num-htvc-threads THREADS \
    --minimum-mapping-quality 1 \
    --min-base-quality-score 35 \
    --f1r2-median-mq 1 \
    --base-quality-score-threshold 35 \
    --allow-non-unique-kmers-in-ref \
    --enable-dynamic-read-disqualification-for-genotyping \
    --recover-all-dangling-branches \
    --genotype-germline-sites \
    --initial-tumor-lod -5 \
    --tumor-lod-to-emit -5 \
    --pruning-lod-threshold -9 \
    --active-probability-threshold 0.00005 \
    --max-mnp-distance 0 \
    --mutectcaller-options '-A AssemblyComplexity -min-dangling-branch-length 0 -max-reads-per-alignment-start 0' \
    --filter-reads-too-long

samtools index PARABRICKS_OUTPUT.bam
bcftools index --tbi PARABRICKS_OUTPUT.vcf.gz
```

The GATK command uses `AssemblyComplexity` and `TandemRepeat` annotations for
SVC model features.
The current NVIDIA Parabricks Mutectcaller documentation lists
`AssemblyComplexity`, but does not list `TandemRepeat`, among the annotations
accepted through `--mutectcaller-options`.
Validate the caller outputs before treating their annotation sets as
equivalent.

See the [NVIDIA Parabricks Mutectcaller user guide](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/mutectcaller)
for other parameters.

Use the resulting VCF and realigned BAM with
`filter_variants tumor-only-wgs`, together with the indexed PON VCF generated
as described in [Generate a Panel of Normals](usage-generate-pon.md).
