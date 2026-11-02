# Troubleshooting

## SID pool size limit

**Symptom:** A sample sheet with more than 256 samples cannot be processed in a single Demux run, regardless of whether indexing is single or dual.

**Cause:** Demux 1.2.0 uses 8-bit identifiers for output routing, which limits a run to 256 sample entries.

**Resolution:** Use a sample sheet with 256 or fewer samples for Demux 1.2.0.
For more than 256 samples, split the sample sheet into batches of 256 or fewer and rerun Demux on the `raw_failed` reads from each preceding run.
This workaround can bias assignments toward earlier batches and should be used only until larger sample sets are supported.
Support for more than 256 samples per run is planned.

## YS non-specific SID assignments

This existing YS behavior was identified during Demux 1.2.0 validation.
It is not related to the custom adapter support added in 1.2.0.

**Symptom:** YS can assign a read to a sample when it finds a non-specific SID match, even if the read has no complete adapter structure.
This can appear as an unexpectedly high number of partial reads or as an elevated `sid_discordant_reads` / `sid_discordant_reads_percentage` rate when non-specific matches at the two ends disagree.
For context, a 65 000-read control containing random A/C/G/T sequences, with read lengths sampled from representative duplex and simplex test inputs rather than a single fixed length and a 105-SID pool, produced 320 YS partial assignments (0.49%).
If those assignments were evenly distributed across the 105 SIDs, the average would be approximately 0.0047% of input reads per SID; this is an average, not a separately measured per-SID rate.
Because the control sequences have no adapter structure, every assignment is a false positive by construction.
The control measured partial assignments only; a separate full-read false-assignment rate was not measured. Full-read false assignments require non-specific matches at both ends and should be much rarer. Under an independence assumption, a two-end rate scales approximately with the square of a per-end probability, but the partial-assignment rate above is a whole-read rate and cannot be squared directly.

**Cause:** A single SID match can be sufficient for YS assignment, even without a corroborating SID at the other end.
Larger SID pools increase non-specific matches across designs.
YS is comparatively less specific because of both its adapter structure and matching algorithm, so the effect is more noticeable.
As with other adapter designs, use only the SIDs of interest where possible.

**Resolution:** If YS produces unexpectedly high rates of assignment, partial reads, or discordant SIDs, use `--read-length-mode full-only` to discard them.
Otherwise, review partial-read output and the discordant-SID metrics before downstream analysis.

## SBX-D long-read trimming

This limitation has been present since the first public Demux release and was identified during Demux 1.2.0 validation.

**Symptom:** The run-level `metrics/read_length_distr/total_trimmed_duplex_passing_read_length_distr.tsv`
(or the corresponding `sample_trimmed_duplex_passing_read_length_distr.tsv` file) can show a small inflection around 427 to 430 bases in the trimmed, passing-read distribution.
This is the characteristic signature of the issue associated with the heuristic's approximately 450-base threshold; it is not necessarily a discontinuity exactly at 450 bases.
Long consensus reads can be over-trimmed when the double-read heuristic is used.

**Cause:** When a consensus is longer than 450 bases, a legacy double-read heuristic searches backwards for a start-adapter match.
An unrelated or random near-match in a legitimate long consensus can be selected as the trim boundary.

**Resolution:** No configuration workaround is available in Demux 1.2.0.
Review read-length distributions and confirm trim boundaries for affected long reads before downstream analysis.

## Custom adapter design limitations

Custom adapter support is new in Demux 1.2.0.
These limitations apply to user-supplied designs in *adapter-config.json5*.
These custom-configuration cases do not affect bundled adapter designs.

### Short fixed terminal sequences

**Symptom:** Demux exits during design loading when a custom duplex adapter has a fixed start or end sequence that is 6 or 7 bases long.
Short terminal sequences can also produce non-specific trimming.

**Cause:** The current matcher evaluates only the shared untrimmed R1 and R2 sequences.
It does not use the sequence immediately before the start adapter or after the end adapter to confirm a terminal match.
The design loader rejects 6- and 7-base terminals, and the fixed edit-distance tolerance makes short terminal matches susceptible to random near-matches.

**Resolution:** Do not use fixed terminal sequences shorter than 8 bases.
Sequences of at least 8 bases can load, but successful loading does not establish adequate specificity.
Use a bundled design when possible, or validate a custom design with representative positive and high-diversity negative-control input before analysis.

### Duplex-family linker and stem structures

**Symptom:** A custom `duplex`-family design with linker or stem geometry different from the bundled structures produces trimmed reads that are too short or too long without reporting an error.

**Cause:** Demux 1.2.0 applies a hardcoded one-base offset in shared duplex trimming calculations.
The `duplex` DSL explicitly represents the standard one-base linker as `A` after `[ins+]` and its reverse-complement `T` before `[ins-]`.
The hardcoded offset can nevertheless mis-handle custom `duplex_stem` and `duplex_umi` structures with different stem or linker geometry.
Designs with no linker base cannot be represented correctly by the current `duplex` trimming path.

**Resolution:** Use `duplex` with its explicit one-base `A`/`T` linker literals; do not omit them, change that structure, or add a custom spacer length.
Use `duplex_stem` for a longer fixed stem, or `duplex_umi` when the structure includes UMIs, and validate output against representative reads.
Do not use a linkerless `duplex` design in Demux 1.2.0.
Use bundled structures until trimming supports structure-derived offsets.

### Highly complex custom simplex fixed regions

**Symptom:** An exceptionally complex custom `simplex` design can fail during design loading if no usable extension DFA can be built.
If an insert-facing stem can be built only from a shortened sequence, part of that stem can remain attached to the demuxed insert.

**Cause:** Generic simplex matching extends SID matches with bounded DFA classifiers.
Read-facing anchors are limited to 15 bases and retried at shorter lengths when the DFA state space is too large.
Insert-facing stems are attempted at full length and use the same retry, so a shorter DFA can leave the distal stem bases untrimmed.
Bitap is used for end-anchored partial-stem recovery, not as a general replacement for complete-read extension.

**Resolution:** This limitation is unlikely to affect bundled designs or most custom designs.
Validate an exceptional custom design with representative reads and inspect trimmed output.
If design loading fails or insert-facing fixed bases remain, simplify the fixed regions or use a bundled design.

## Very few reads assigned with a custom sample sheet

**Symptom:** You supply your own sample sheet with `--sample-sheet` and the
adapter design is correct, but almost all reads land in `unassigned` and no error
is raised.

**Cause:** One or more `sample_sid` values do not match a canonical HPSID in
either orientation, or the sample sheet mixes SID representations. A legacy SID
that includes adjacent adapter bases is recognized only when those bases match
the selected adapter design's DSL and every SID uses the same orientation and
flank placement.

**Resolution:** Use canonical HPSID sequences for every adapter design (see
[SID orientation](file-formats.md#sid-orientation)).
Demux automatically converts a sample sheet when every SID is a recognized
reverse complement of a canonical HPSID, or when every SID uses the same
DSL-valid legacy representation containing adjacent structural bases. Demux
preserves mixed representations, invalid structural flanks, and unrecognized
custom SID pools as supplied and logs a warning.

## High memory usage or crashes on large inputs

**Symptom:** `demux` consumes more memory than available, or is killed by the OS.

**Cause:** Peak memory scales with `--threads` and `--batch-size`; strand
detection additionally holds the strand index in memory.

**Resolution:** Lower `--batch-size`, reduce `--threads`, or raise
`--workers-per-reads-file` to process fewer input files concurrently. With
defaults, allocate roughly 1 GB per 10 threads. For `demux_strand_index`, building
the human genome index needs ~32 GiB (see
[Performance and cost notes](appendix.md#performance-and-cost-notes)).

## demux_strand_index fails to read the reference FASTA

**Symptom:** `demux_strand_index` errors out reading the `--reference` FASTA.

**Cause:** The reference FASTA is missing its `.fai` index.

**Resolution:** Create the index first (e.g., `samtools faidx reference.fa`) so a
companion `reference.fa.fai` sits next to the FASTA, then re-run.
