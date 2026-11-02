# Building the strand index

`demux_strand_index` builds a strand index — an interleaved Bloom filter of
reference genome k-mers — used for optional, alignment-free
[genome strand detection](appendix.md#strand-detection) of duplex reads in
[`demux`](usage-duplex.md). Input: **reference FASTA** → Output: **`strand.bloom`
index file**.

## How it works

The reference genome is decomposed into k-mers. For each k-mer, its presence is
recorded in a forward-strand Bloom filter and its reverse complement in a
reverse-strand filter. The two filters are **interleaved** so that both strands
for a given k-mer sit adjacent in memory, which speeds up lookups when `demux`
queries the index. During a duplex run, `demux` queries this index to decide
which genomic strand a read originated from, enabling strand-aware indel
left-alignment without changing the consensus sequence.

Expand the details below for the strand-detection model and the main sizing
heuristics.

<details>

<summary>Design details</summary>

Genome strand detection is optional and does not change consensus base calls — it
only affects how indels are left-aligned (and therefore where they are placed).
See [Genome strand detection](appendix.md#strand-detection) for the full
background.

Filter sizes are powers of 2, and entropy-optimal false-positive rates are
`1/2^n`. The default rate `0.015625` (n=6) sits near optimal entropy: for the
human genome it produces a ~2 GiB filter.

</details>

## Example command

```bash
demux_strand_index \
  --reference /path/to/reference_genome.fa \
  --threads 16
```

Creates `strand.bloom` in the current directory. Pass it to `demux` via
`--strand-index`:

```bash
demux \
  --adapter-design-bundle /resources/adapter-design-bundle.zip \
  --adapter-design-name SBX-D \
  --threads 16 \
  --reads /path/to/input_directory \
  --output-dir /path/to/output_directory \
  --strand-index strand.bloom
```

## Input

| Input | Required | Notes |
|:------|:---------|:------|
| Reference FASTA | Yes | Passed with `--reference`. Reference genome FASTA. Must have a companion `.fa.fai` index. |

## Output

| Output | Description |
|:-------|:------------|
| *strand.bloom* | Interleaved Bloom filter written to the current directory (override with `--output-strand-index`). Pass to `demux --strand-index`. |

## Tips for analysis

Expand the note below for guidance on building and tuning the index.

<details>

<summary>Building and tuning the index</summary>

1. Build once per reference genome; reuse across runs.
2. Larger `--kmer-size` improves precision but reduces tolerance to sequencing
   errors.
3. Lower `--max-false-positive-rate` reduces false positives but increases filter
   size.
4. Entropy-optimal FPRs are `1/2^n`. The default `0.015625` (n=6) works well —
   only change it if memory is constrained.
5. Filter sizes are powers of 2. For the human genome the default FPR produces a
   2 GiB filter near optimal entropy.

</details>

Performance and cost notes for simplex, duplex, and strand-index workflows are
described in the [Appendix](appendix.md#performance-and-cost-notes).

## CLI options

This command also accepts the common `-h, --help`, `-v, --version`, and
`-l, --log-level` options.

Parameters in **bold** are required.

### Input options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| **--reference** | Reference FASTA file. Must include a fai index. | A string file path |

### Output options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --output-strand-index | Output filename for the strand index. | A string file path [default `strand.bloom`] |
| --kmer-size | k-mer size for the output strand index. | An integer > 0 and < 32 [default `19`] |
| --max-false-positive-rate | Maximum false positive rate for the output interleaved bloom filter used as index. | A float in range (0, 1) [default `0.015625`] |

### Performance options

| Parameter | Description | Value(s) |
|:----------|:------------|:---------|
| --threads | Number of threads used (0=available hardware threads). | An integer >= 0 [default `1`] |
