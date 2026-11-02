# Usage

This page helps you choose the command guide that matches your library type and
workflow.
`demux` is one command with two configuration-specific usage paths.
Choose the page that matches the adapter design selected by
`--adapter-design-name`.

## Simplex demultiplexing

Demultiplexes, trims, and filters raw _simplex_ sequencing reads.
Selected when `--adapter-design-name` names a simplex design (YS, YSU,
SIMPLEX-10X).
Input: **raw FASTQ/RDB** → Output: **per-sample FASTQ**.

```bash
demux --adapter-design-name YSU --reads raw_reads/ --output-dir out/
```

To run simplex demultiplexing, go to
[Simplex demultiplexing](usage-simplex.md).

## Duplex demultiplexing

Demultiplexes, trims, filters, and computes intramolecular consensus for raw
_duplex_ sequencing reads.
Selected when `--adapter-design-name` names a duplex design (SBX-D, SBX-DM,
SBX-FAST).
Input: **raw FASTQ/RDB** → Output: **per-sample consensus FASTQ**.

```bash
demux --adapter-design-name SBX-D --reads raw_reads/ --output-dir out/
```

To run duplex demultiplexing, go to
[Duplex demultiplexing](usage-duplex.md).

## Building the strand index

Builds an interleaved Bloom filter of reference genome k-mers for optional,
alignment-free [genome strand detection](appendix.md#strand-detection)
of duplex reads.
Input: **reference FASTA** → Output: **`strand.bloom`**.

```bash
demux_strand_index --reference reference.fa
```

To build the strand index, go to
[Building the strand index](usage-strand-index.md).
