# Troubleshooting

Common issues when running the STR Caller and how to resolve them.
Issues are ordered by how frequently users hit them.

## A locus is flagged LowDepth or LowConfidence

**Symptom:** A locus is filtered `LowDepth`, or a call is filtered `LowConfidence`.

**Cause:** Too few reads were available at the locus for a reliable call. This is common when `--skip-r1-r2-decoding` is set, which reduces the number of usable reads per site. If the site is `LowConfidence` and the variant type is `RepeatFXN`, this maybe a case where the fraction of spanning reads is outside the range to make a confident call.

**Resolution:** Inspect the JSON output for the per-locus repeat counts and the distribution of spanning/flanking/in-repeat reads to judge the evidence. Avoid `--skip-r1-r2-decoding` unless required, and confirm the locus has adequate coverage in the input BAM. For `RepeatFXN` variant type, this means we do not have suffcient read evidence to make a confident call.

## A call looks wrong and you need to inspect the evidence

**Symptom:** A genotype call is unexpected and the VCF alone does not explain it.

**Cause:** The VCF summarizes the call but not the underlying read support and realignments.

**Resolution:** Use the JSON output for repeat counts and the spanning/flanking/in-repeat read distribution. Set `--output-bam` to emit the recruited reads and their realignments; without `--write-reads-as-unmapped`, the BAM retains the original reference coordinates of the recruited reads for browsing.

## The wrong aligner is assumed for the BAM

**Symptom:** Off-target read recruitment behaves unexpectedly for a Giraffe-aligned BAM.

**Cause:** When `--aligner` is omitted, the aligner is inferred from the BAM `@PG` header, and off-target recruitment differs between `bwa` and `giraffe`. Detection can be incorrect for some BAMs.

**Resolution:** Set `--aligner` explicitly (`bwa` or `giraffe`) to match how the BAM was produced.

## The run stops because output files already exist

**Symptom:** The caller exits rather than overwriting existing output.

**Cause:** By default the caller does not overwrite existing output files.

**Resolution:** Pass `--overwrite`, or choose a fresh `--output-dir`.
