# Appendix

## Concepts and terminology

| Term | Definition                                                                                                                                                                                                                                            |
|:-----|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Short Tandem Repeat (STR) | A short sequence of DNA bases repeated in tandem at a specific location in the genome. STRs typically have 2-6 base-pair unit lengths and can vary in the number of repeat units among individuals.                                                   |
| Spanning reads | Reads that overlap both flanking regions of a repeat region, providing the strongest support for estimating an exact repeat count.                                                                                                                    |
| Flanking reads | Reads that overlap only one of the two flanking regions of a repeat region.                                                                                                                                                                           |
| In-repeat reads | Reads that cover only the repeat region without any flanking region overlap. Generally observed when repeat lengths expand beyond the read length, these reads can only be used to improve confidence interval estimates and not exact repeat counts. |
| Local realignment blueprint | A representation of the structure of a repeat region, including the arrangement of repeat units and flanking sequences. Blueprints are used to construct local realignments for aligning reads.                                                       |
| VariantType | The classification of a variant as either a common repeat (`Repeat`) or a rare repeat (`RareRepeat`). This affects how reads are recruited and analyzed.                                                                                              |
| Off-target regions | Genomic regions where reads from expanded repeats can misalign. These are used to recruit additional reads for rare-repeat analysis.                                                                                                                  |
| Weighted purity score | A metric evaluating what fraction of the read sequence consists entirely of the repeat unit, based on the proportion of matching bases weighted by base-call quality scores.                                                                          |
| Genotype confidence interval | A range of repeat sizes considered plausible for a genotype, given the observed read data. It measures uncertainty in the genotype estimate.                                                                                                          |

## Experimental features and limitations

The pileup visualization, projected-read IGV files, and STR-region variant VCF are experimental outputs that have undergone limited testing.
Use them to inspect supporting evidence and sequence variation, and interpret results cautiously for loci or data types outside the evaluated examples.

The STR-region variant VCF reports SNVs and indels inside and around an STR locus. For this experimental feature, variant calling relies solely on reads anchored by flanking sequence. Variants present entirely within in-repeat reads will not be called.
The `LocusSequenceVariant` field in the repeat-count VCF is a locus-level indicator that points to the individual sequence-level records in the STR-region variant VCF.

The caller supports BAMs produced by BWA-family aligners (for example, BWA-MEM or BWA-MEM2) and `vg giraffe`.
Read interpretation and off-target recruitment can differ between aligners, so use the aligner that produced the BAM or set `--aligner` explicitly to the exact supported value that matches your BAM when auto-detection is not reliable.

## Performance and cost notes

Catalog entries are processed in parallel, so runtime scales with the number of loci and their read depth, and decreases roughly linearly with `--threads`.
`--region-padding` widens the recruitment window per locus, trading additional read processing for sensitivity to reads anchored further from the repeat.

Running with a 30X BAM (~700M reads) with 40 cataloged regions, the following performance was observed, estimating 200 GB storage and an AWS on demand instance of type m5.4xlarge (16 vCPU, 64 GiB memory) at $0.768 per hour down to m5.large (2 vCPU, 8 GiB memory) at $0.096 per hour:

| Num Threads | Time (s) | Cost $(USD) | Instance Type |
|:-------------|:----------|:-------------|:---------------|
| 1           | 39       | 0.06        | m5.large      |
| 2           | 20       | 0.03        | m5.large      |
| 4           | 11       | 0.04        | m5.xlarge     |
| 8           | 6        | 0.04        | m5.2xlarge    |
| 16          | 4        | 0.05        | m5.4xlarge    |
