# Appendix

## Concepts and terminology

| Term | Definition |
|:---|:---|
| Duplex consensus | Duplex sequencing reads both strands of a DNA fragment, producing two related reads (R1 and R2). A consensus base is derived at each position by comparing the two strands. |
| Discordant base | A position where the R1 and R2 bases disagree. When R1 and R2 disagree, the R1 base is chosen by default during consensus. |
| Pangenome (graph) alignment | An alignment of a read against a reference pangenome graph rather than a single linear reference. A graph can represent known variation as alternate paths, so a read carrying that variation aligns without mismatches. Input BAM files must be aligned with [Giraffe](https://github.com/vgteam/vg) using the `--add-graph-aln` parameter, which records the graph alignment in the `GR` tag. |
| GR tag | A per-read SAM tag written by Giraffe that encodes the read's alignment to the pangenome graph as a *difference string* (`cs`-style). The Pan-Genome Consensus Caller translates the `GR` tag into the graph base aligned at each read position. A `GR` value of `*` indicates no graph alignment, and the read is passed through unchanged. |
| YC tag | A per-read SAM tag describing the duplex structure of the read, including which positions are simplex, which are duplex matches, and which are discordant. For a full description of the `YC` tag, see [Interpreting SBX Duplex (SBX-D) Reads](https://xoos.gitbook.io/xoos). |
| Difference string | A compact string encoding of an alignment as a sequence of operations (matches, substitutions, insertions, and deletions) relative to the aligned sequence. The `GR` tag stores the graph alignment in this form. |

## Performance and cost notes

Processing scales with the number of threads set by `-t`.
When more than one thread is used, output records are emitted in a non-deterministic order and the output BAM is unsorted; sort it with `samtools sort` before indexing.
Running single-threaded preserves input order.
