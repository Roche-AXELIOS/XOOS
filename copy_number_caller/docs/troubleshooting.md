# Troubleshooting

Common issues and their resolutions, ordered by how often users hit them.

## Empty output files from low-coverage input

**Symptom:** Copy Number Caller exits with code 0 but writes structurally valid, data-empty output files (no segments or variants).

**Cause:** The input BAM has a median coverage below 10x, or, in somatic modes, the VCF has too few germline SNPs or zero variants after filtering. The module skips calling to avoid numerical errors during GC correction.

**Resolution:** This is expected behavior for insufficient-coverage input. Confirm the input BAM meets the coverage requirement and, for somatic modes, that the VCF contains enough germline heterozygous SNPs. Re-run with adequate input.

## Missing common duplications from pangenome-aligned BAMs

**Symptom:** Germline CNV calls from a giraffe pangenome-aligned BAM miss common duplications.

**Cause:** Reads that align to the duplicated copy in the pangenome graph can fail to project to the copy in the linear reference, so the duplicated region appears at normal copy number. BWA-aligned BAMs are not affected.

**Resolution:** Run giraffe with `--max-multimaps` so the rescue and tag method recovers these reads. This is enabled automatically for CPU-based vg giraffe in the XOOS pipeline. See the [`germline-wgs` usage page](usage-germline-wgs.md#how-it-works) for details.

## Blank output files on Apple M-series machines

**Symptom:** Blank output files are produced when running on an Apple M-series (M-chip) machine, both when running the binary directly and through Docker.

**Cause:** The third-party scientific computing libraries used during GC correction do not produce correct results on this architecture.

**Resolution:** Run Copy Number Caller on a supported x86_64 host (see the recommended system requirements in the module README).

## Preventing oversegmentation

**Symptom:** The output contains excessive numbers of small segments, particularly in noisy samples such as FFPE.

**Cause:** High log-ratio noise can create spurious segment boundaries and fragment the genome.

**Resolution:** Use `--segment-undoing` to merge adjacent segments with similar log-ratios, adjusting its strategy and factors as needed.

Applies to: [`tumor-normal-wgs`](usage-tumor-normal-wgs.md) and [`tumor-only-wgs`](usage-tumor-only-wgs.md).

For details about segment undoing strategies and tuning, see [Oversegmentation](appendix.md#oversegmentation).

## Predicting SCNV in low SCNV signal tumors

**Symptom:** Copy-number calls or purity/ploidy estimates are unreliable, often alongside weak log-ratio and BAF signals.

**Cause:** Low tumor purity and/or low SCNV burden weakens the signal used to estimate purity, ploidy, and copy number.

**Resolution:** Treat calls from tumors with predicted purity below 30% as generally unreliable. Inspect the log-ratio and BAF channels regardless of predicted purity; for suspected high-purity tumors with low SCNV burden, provide a somatic flagged VCF to estimate purity from somatic variants.

Applies to: [`tumor-normal-wgs`](usage-tumor-normal-wgs.md) and [`tumor-only-wgs`](usage-tumor-only-wgs.md).

For background on low SCNV signal, its effect on purity/ploidy estimates, and the tumor-normal somatic-variant fallback, see [Low SCNV signal tumors](appendix.md#low-scnv-signal-tumors).

## Purity and ploidy ambiguity in whole-genome doubled tumors

**Symptom:** The estimated purity/ploidy solution or resulting calls may conflict with expectations for a genome-doubled tumor.

**Cause:** Multiple purity/ploidy combinations can fit the data similarly, and the diploid prior favors lower-ploidy solutions.

**Resolution:** Inspect the log-ratio and BAF channels and the purity/ploidy likelihood grid. If needed, use `predict-somatic-cnv` with an alternative solution to regenerate calls.

Applies to: [`tumor-normal-wgs`](usage-tumor-normal-wgs.md) and [`tumor-only-wgs`](usage-tumor-only-wgs.md).

For background on purity/ploidy ambiguity and whole-genome-doubling checks, see [Whole genome doubling](appendix.md#whole-genome-doubling).
