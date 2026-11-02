# File formats

This page is the reference for the Demux-specific input and output formats: the
sample sheet, the adapter design bundle, output read-name annotations, and the
output directory layout. For how these formats are used in practice, see
[Simplex demultiplexing](usage-simplex.md#output),
[Duplex demultiplexing](usage-duplex.md#output), and
[Building the strand index](usage-strand-index.md#output). Metric-file schemas
live on the [Metrics](metrics.md) page. Standard containers (FASTQ, FASTA) used
in the ordinary way are not re-documented here.

## Sample sheet format

Applies to: `demux` simplex, `demux` duplex.

The sample sheet is a CSV with two columns, `sample_name` and `sample_sid`:

```text
sample_name,sample_sid
HG001,CTATTGAAGGCT
HG002,GGAGTAACATCC
```

Supply a sample sheet with `--sample-sheet` to restrict SID lookup to these
samples. If omitted, the adapter bundle's default SID pool is used, which can
include samples not associated with the run.

### Dual-SID designs

For dual-SID simplex designs such as `I2X-S-DX`, write `sample_sid` as
`sid_5p-sid_3p` in the orientation of the corresponding DSL placeholders:

```text
sample_name,sample_sid
sampleA,AGCCTTCAATAG-GGATTATGGCAA
sampleB,CTATATGAGGCG-GGATTATGGCAA
```

Same-SID rows without a hyphen are normalized internally as `sid-sid`. `I2X-S-DX`
requires `--sample-sheet`; bundled SID pair definitions are not provided.

### SID orientation

{% hint style="info" %}
Enter the canonical SID sequence for every adapter design.
The adapter design applies the required orientation during matching.
{% endhint %}

Demux uses [*sids/hpsid.fa*](../resources/adapter-design-bundle/sids/hpsid.fa)
as the single canonical list of SIDs for the built-in adapter designs.
Every sequence in this FASTA is stored in canonical orientation:

```text
>HPSID_001_L1
CTATTGAAGGCT
```

Use the canonical sequence in the sample sheet, regardless of the order or
orientation in which the SID occurs in a read:

```text
sample_name,sample_sid
HG001,CTATTGAAGGCT
```

The adapter design's domain-specific language (DSL) controls the sequence that
Demux matches at each SID position.
For this example, `[sid:sid+]` matches the canonical `CTATTGAAGGCT`, and
`[sid:sid-]` matches its reverse complement, `AGCCTTCAATAG`.

#### YS example

The built-in YS design declares the reverse-complemented SID first and the
canonical SID second:

```text
CAACAAGAGTCTTTT[sid:sid-]GCGT[ins]ACGC[sid:sid+]TTTGTCGTGTAGGTGAAGAGC
```

The corresponding SID regions in a read containing `HPSID_001_L1` are:

```text
...AGCCTTCAATAGGCGT[insert]ACGCCTATTGAAGGCT...
   [sid:sid-]                   [sid:sid+]
```

#### SBX-D example

The built-in SBX-D design declares the canonical SID first and the
reverse-complemented SID second:

```text
CAACAAT(GCACCTGACGTACGCT[ins+]A[sid:sid+])ATGTATG([sid:sid-]T[ins-]AGCGTACGTCAGGTGC)TGAAGAGC
```

The corresponding SID regions in a read containing `HPSID_001_L1` are:

```text
...[insert+]ACTATTGAAGGCT...ATGTATG...AGCCTTCAATAGT[insert-]...
             [sid:sid+]                 [sid:sid-]
```

For either design, Demux assigns the example read to `HG001`, writes it under
the *HG001/* output directory, and reports `CTATTGAAGGCT` in the
`index_sequence` row of `sample_stats.tsv`.
Without a sample sheet, Demux uses the canonical FASTA record name instead, so
the same SID is assigned to `HPSID_001_L1` and written under the
*HPSID_001_L1/* output directory.

{% hint style="info" %}
For compatibility, Demux normalizes a sample sheet when every SID is either the
reverse complement of a canonical HPSID or the same legacy representation of a
canonical HPSID with adjacent structural bases included. Legacy structural
bases must match the selected adapter design's DSL, including their orientation
and placement. Demux converts the complete pool to canonical orientation and
logs a warning. It preserves mixed representations, invalid structural flanks,
and unrecognized custom pools as supplied and logs the applicable warning.
Automatic conversion is for legacy compatibility only. To remove its warning,
replace the converted sample-sheet values with their canonical 12-base SIDs
from *sids/hpsid.fa*.
{% endhint %}

### Advanced: custom SID orientation

Normalization is driven by the SID FASTA declared by the selected adapter
design; there is no separate normalization flag. If a custom adapter intentionally
uses a different SID orientation, define that orientation as the canonical one in
its own SID FASTA and use the matching DSL orientation. For example, this design
swaps the two YS SID orientations:

```json5
{
  "adapter_designs": {
    "YS-CUSTOM": {
      "type": "simplex",
      "dsl": "CAACAAGAGTCTTTT[sid:sid+]GCGT[ins]ACGC[sid:sid-]TTTGTCGTGTAGGTGAAGAGC",
      "adapters": {
        "sid": { "sequence": "sids/custom-sids.fa" }
      }
    }
  }
}
```

In this example, `sids/custom-sids.fa` must contain the SID sequences in the
orientation used by the sample sheet (for example, `AGCCTTCAATAG` rather than
the built-in canonical `CTATTGAAGGCT`):

```text
>HPSID_001_L1
AGCCTTCAATAG
```

Because those values are canonical for `YS-CUSTOM`, Demux leaves them unchanged
and applies the custom DSL. This is an advanced escape hatch for intentionally
custom adapter work; canonical sample-sheet values remain the recommended form
for built-in designs.

Back to [Simplex demultiplexing input](usage-simplex.md#input) and
[Duplex demultiplexing input](usage-duplex.md#input).

## Adapter bundle structure

Applies to: `demux` simplex, `demux` duplex.

The adapter design bundle (`--adapter-design-bundle`) is a folder or ZIP
containing FASTA files and an `adapter-config.json5` that describes the
available designs:

```text
adapter_design_bundle/
├── adapter_name/
│   ├── adapter_feature.fa
└── adapter-config.json5
```

Each design is expressed as a short read-structure string (a domain-specific
language, or DSL) plus a small set of per-barcode settings.
You describe the adapter layout once, and Demux derives the individual lookup
tables and search windows from it.
This keeps a design readable end to end and lets you add a new adapter by
editing config alone, without a code change.

### Custom adapter design with DSL specification

Expand the details below for the DSL authoring model and how Demux compiles it
into the runtime representation.

<details>

<summary>DSL concepts and compilation model</summary>

The adapter design DSL lets you describe an adapter as a single string, written
5′→3′ the way it appears on a spec sheet.
The structure *is* the string: element order, orientation, and linkage are read
directly from the text, so a design can be authored, shared, and compiled
without hand-assembling parallel arrays.

**Two forms, one model:**

- **Authored** — what you write: a `dsl` string plus an `adapters` map (id →
  FASTA path and settings).
- **Compiled** — what Demux consumes: a `read_structure` token array plus the
  `adapters` map with each entry's `type` made explicit.

```text
DSL (you write) -> { read_structure, adapters } (Demux reads)
```

A design file accepts two top-level shapes: the preferred wrapper, which maps
design names under `adapter_designs`, or a flat map of design name to design
object.
An optional `syntax` field on an authored design selects the authoring grammar
and defaults to the DSL.
The single-design objects shown below use the flat-map shape for brevity.

**DSL syntax:**

- **Fixed sequence** — literal bases, written inline in uppercase (for example
  `CAACAAT`).
- **Special elements** in square brackets: `[sid:<id>]`, `[umi:<id>]`, and
  `[ins]`.
  The first three characters select the kind: `sid` (sample index, drives
  demultiplexing), `umi` (unique molecular index), and `ins` (the biological
  insert, a reserved keyword with no id).
- **IDs** follow `^[a-z][a-z0-9_]+$` and reference a FASTA of candidate
  barcodes.
  Elements sharing a base id (for example `sid+` and `sid-`) are linked and
  must match by index for a valid result; distinct ids (for example `umi1`,
  `umi2`) are independent.
- **Orientation** — a trailing `+` / `-` declares strand directly.
  `sid-` means reverse-complemented relative to its linked FASTA.
  `+` and no marker are equivalent; `+` is recommended for clarity.
- **Duplex fold** — parentheses `(...)` mark each arm of a duplex hairpin.
  The two arms must be reverse complements of each other, and the fixed
  sequence between them is the loop.

**Compiled conventions:**

The `read_structure` token list is self-describing with no per-token metadata.
Token case is the class (`^[ATCG]+$` for a literal sequence;
`^[a-z][a-z0-9_]+$` for an id or the reserved `insert` keyword), the optional
`+` / `-` suffix carries orientation, matching base names indicate linkage, and
each duplex arm becomes a nested array.
The `adapters` map is keyed by id base name and each entry declares a `type` of
`sid`, `umi`, or `fixed`.
`sid` / `umi` entries supply the FASTA `sequence` and any search parameters; a
`fixed` entry is optional and only attaches per-marker search settings to a bare
fixed literal that already appears in `read_structure`.
The `insert` keyword is read-structure-only and never has an `adapters` entry.

**Example — SBX-D (duplex):**

Authored:

```json
{
  "SBX-D": {
    "type": "duplex",
    "dsl": "CAACAAT(GCACCTGACGTACGCT[ins+]A[sid:sid+])ATGTATG([sid:sid-]T[ins-]AGCGTACGTCAGGTGC)TGAAGAGC",
    "adapters": { "sid": { "sequence": "sids/hpsid.fa" } }
  }
}
```

Compiled:

```json
{
  "SBX-D": {
    "type": "duplex",
    "read_structure": [
      "CAACAAT",
      ["GCACCTGACGTACGCT", "insert+", "A", "sid+"],
      "ATGTATG",
      ["sid-", "T", "insert-", "AGCGTACGTCAGGTGC"],
      "TGAAGAGC"
    ],
    "adapters": { "sid": { "type": "sid", "sequence": "sids/hpsid.fa" } }
  }
}
```

**Compilation steps:**

1. Tokenize into bracketed elements and the bare uppercase runs between them.
2. Classify each bracket by its 3-character prefix.
   `ins` is the reserved keyword (`[ins]`, `[ins+]`, `[ins-]`; no id, no
   colon), while `sid` / `umi` require `:<id>` with an optional trailing `+` /
   `-`.
3. Emit `read_structure`.
   Fixed sequences pass through, `sid` / `umi` ids become bare `id±` tokens,
   and `ins` becomes `insert`.
4. Build the `adapters` skeleton keyed by id base name with the DSL-derived
   `type` (`sid` / `umi`).
5. Overlay the authored `adapters` block onto that skeleton, supplying each
   id's FASTA `sequence` and search parameters, and permitting extra
   `type: fixed` entries only when they anchor a fixed literal present in
   `read_structure`.
6. Validate.
   Every DSL id must have a matching authored entry (and vice versa), the
   declared and DSL types must agree, `insert` can never be an id, arm symmetry
   must hold for duplex designs.

**Adapter entry fields:**

Each `sid` / `umi` entry carries the FASTA `sequence` and, at minimum, its
`type`.
It can also carry these optional search parameters: `max_edit_distance`, the
strand-specific wiggle windows `wiggle+` / `wiggle-`, the extension lengths
`additional_prefix_length` / `additional_suffix_length`, and the bait lengths
`bait_length+` / `bait_length-`.
A `fixed` entry can carry a single `wiggle` window (`left` / `right`) to tune
the search around a fixed literal.

{% hint style="info" %}
`max_edit_distance`, `wiggle` (`wiggle+` / `wiggle-` and the fixed-entry
`wiggle`), `bait_length+` / `bait_length-`, and
`additional_prefix_length` / `additional_suffix_length` are legacy search knobs
that are being phased out in favor of inferential methods.
Once inference fully replaces them, these fields will be removed.
The `type` and the FASTA `sequence` are intrinsic to the molecule and remain
required.
{% endhint %}

</details>

### Adapter config file structure

Expand the details below for the supported top-level fields and per-design
configuration fields.

<details>

<summary>adapter-config.json5 fields</summary>

**Top-level fields:**

- `adapter_designs`: Map of design name to design object (for example `"SBX-D"`,
  `"YSU"`).

**Adapter design fields:**

- `type`: Adapter family.

| Type | Structure | Notes |
|:-----|:----------|:------|
| `duplex` | `[N1]([N2+][ins+]A[sid:sid+])[N3]([sid:sid-]T[ins-][N2-])[N4]` | `()` arms are reverse complements. |
| `duplex_stem` | `[N1]([N2+][ins+][N3+][sid:sid+])[N4]([sid:sid-][N3-][ins-][N2-])` | Adds fixed stem sequences `[N3]`. |
| `duplex_umi` | `[N1]([N2+][umi:umi1+][ins+][umi:umi2-][N3+][sid:sid+])[N4]([sid:sid-][N3-][umi:umi2+][ins-][umi:umi1-][N2-])` | Duplex stem with `[umi]` support. |
| `YSU` | `[N1][sid:sid-][N2+][umi:umi1+][N3+][ins][N3-][umi:umi2-][N2-][sid:sid+][N4]` | |
| `YS` | `[N1][sid:sid-][N2+][ins][N2-][sid:sid+][N4]` | |
| `simplex` | `[N1][sid:sid-][N2+][ins][N2-][sid:sid+][N4]` | Same layout as `YS`, but distinct for downstream processing. |
| `simplex-10x` | `[N1][sid:sid-][N2+][ins]` | No adapters on the 3′ end. |

`Nx` is a numbered fixed sequence.

- `dsl`: The read structure, left to right in read order.
  Literal bases are written directly; barcodes are written as
  `[sid:<id>±]` and `[umi:<id>±]`; the sequenced insert is `[ins]`.
  The `+` or `-` suffix declares the orientation of each barcode occurrence;
  Demux does not derive SID orientation from the adapter `type`.
- `adapters`: Map of barcode id (the `<id>` used in the DSL) to its settings.

**Barcode settings** (within `adapters`):

- `sequence`: Path to the FASTA file of candidate barcodes, relative to the
  bundle root.
- Legacy settings that are being phased out: `max_edit_distance`, `wiggle+` /
  `wiggle-`, `bait_length+` / `bait_length-`,
  `additional_prefix_length`, and `additional_suffix_length`.

**Example snippet:**

```json
{
  "adapter_designs": {
    "YS": {
      "type": "YS",
      "dsl": "CAACAAGAGTCTTTT[sid:sid-]GCGT[ins]ACGC[sid:sid+]TTTGTCGTGTAGGTGAAGAGC",
      "adapters": {
        "sid": {
          "sequence": "sids/hpsid.fa",
          "max_edit_distance": 2,
          "wiggle+": { "left": 32, "right": 2 },
          "wiggle-": { "left": 2, "right": 17 }
        }
      }
    }
  }
}
```

To add your own adapter, copy an existing design of the same family, point
`sequence` at your barcode FASTA, adjust the literal bases and search settings,
and give it a new name under `adapter_designs`.
Select it at run time with `--adapter-design-name`.

</details>

Back to [Simplex demultiplexing input](usage-simplex.md#input) and
[Duplex demultiplexing input](usage-duplex.md#input).

## Output read-name formats

Applies to: `demux` simplex, `demux` duplex.

Demux annotates each output read name with the SID and UMI it detected. The exact
format depends on whether the input was RDB or FASTQ.

### Read name format for FASTQ output with RDB input

When input files are in RDB format (a binary format), output FASTQ read names are
regenerated. Four RDB input flavors are supported — Gamma2, Gamma1, DataPath2, and
HTP1 — and are auto-detected by the reader. Read names use this format:

```text
{date}:{sequence}:Q{queue_name}:R{run_num}:{cycle_id}:{b64-UID}:{sidUmiBitflag}|{sidId}{optionalUMI}
```

| Field | Description |
|:------|:------------|
| `date` | Date of processing in YYYYMMDD format. |
| `sequence` | Sequence number of the read in the input file. |
| `queue_name` | Name of the queue the read was processed in. `*` if not specified or identified. |
| `run_num` | Run number of the read. `*` if not specified or identified. |
| `cycle_id` | Batch index from the RDB block header, used to distinguish reads from different blocks. |
| `b64-UID` | Base64-encoded unique identifier of the read. Encoded using `0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+_`, 0-padded for uniform length. Constant size is 6 tokens (6 tokens × 6 bits = 36 bits, sufficient to exceed u32 range; upper bits are always 0). |
| `sidUmiBitflag` | 4-bit bitfield printed as a base-10 integer (0–15) indicating the presence of SIDs and UMIs. Bits: 5′ SID found `0b1000` (8), 3′ SID found `0b0100` (4), 5′ UMI found `0b0010` (2), 3′ UMI found `0b0001` (1). |
| `sidId` | Assigned SID for the read. |
| `optionalUMI` | UMI sequences if present, formatted as `:{umi5p}:{umi3p}`. |

### Read name format for FASTQ output with FASTQ input

When input files are in FASTQ format, Demux preserves the original read names and
appends the same SID/UMI annotations described above (`sidUmiBitflag`, `sidId`,
and `optionalUMI`) to the end of each read name.

Back to [Simplex demultiplexing output](usage-simplex.md#output) and
[Duplex demultiplexing output](usage-duplex.md#output).

## Output directory and filename structure

Applies to: `demux` simplex, `demux` duplex.

The output directory (`--output-dir`) contains one subdirectory per sample, named for
the sample. Each holds that sample's FASTQ files, named
`<sample_name>-<part_index>.fastq[.gz|.zst]`, where `part_index` is a 1-based,
zero-padded integer (6 digits by default) that uniquely identifies each output
file within a sample. Part indices are contiguous, assigned in input file order,
and files from the same input file use contiguous indices.

The output directory also contains a `metrics/` directory with the files described
in [Metrics](metrics.md).

If `--write-failed-reads` is set, a `raw_failed` sample directory holds
unassigned reads that failed before or during SID assignment, named
`raw_failed-<part_index>.fastq[.gz|.zst]`.

Each sample subdirectory contains a `file_provenance.tsv` that maps each output
FASTQ file to the input file it was derived from (columns: `output_filename`,
`input_filepath`). This is useful for troubleshooting.

When an explicit `--sample-sheet` is provided and a sample listed in it receives
no reads across all inputs, it still produces its per-sample directory containing
a single valid but empty FASTQ (`<sample_name>-000001.fastq[.gz|.zst]`). Because
no input contributed reads to that file, its `file_provenance.tsv` row records
`NA` in the `input_filepath` column.

{% hint style="info" %}
**Determinism:** With `--writers-per-sample 1` (the default), output is
fully deterministic — the same inputs produce the same filenames, read order, and
contents. With more than one writing thread, read order within a sample can vary
between runs due to concurrent writes.
{% endhint %}

A truncated example, with a single writing thread, `--write-failed-reads` set,
and two samples `HG0001` and `HG0002`:

```text
output_directory/
├── metrics/
│   ├── run_stats.tsv
│   ├── sample_stats.tsv
│   ├── sample_assignment_metrics.tsv
│   └── read_length_distr/
│       ├── total_read_length_distr.tsv
│       ├── raw_failed_read_length_distr.tsv
│       ├── total_untrimmed_simplex_passing_read_length_distr.tsv
│       ├── total_trimmed_simplex_passing_read_length_distr.tsv
│       ├── total_endadapter_position_distr.tsv
│       └── HG0001_read_length_distr/
│           ├── sample_untrimmed_simplex_passing_read_length_distr.tsv
│           ├── sample_trimmed_simplex_passing_read_length_distr.tsv
│           └── sample_endadapter_position_distr.tsv
├── HG0001/
│   ├── HG0001-000001.fastq.gz
│   ├── HG0001-000002.fastq.gz
│   ├── HG0001-000003.fastq.gz
│   └── file_provenance.tsv
├── HG0002/
│   ├── HG0002-000001.fastq.gz
│   ├── HG0002-000002.fastq.gz
│   └── file_provenance.tsv
└── raw_failed/
    ├── raw_failed-000001.fastq.gz
    ├── raw_failed-000002.fastq.gz
    └── file_provenance.tsv
```

With `--read-length-mode all-split`, full and partial reads are written to
separate subdirectories:

```text
out-dir/
├── <sample>/
│   ├── full/
│   │   ├── <sample>-000001.fastq.gz
│   │   └── file_provenance.tsv
│   └── partial/
│       ├── <sample>-000001.fastq.gz
│       └── file_provenance.tsv
```

Back to [Simplex demultiplexing output](usage-simplex.md#output) and
[Duplex demultiplexing output](usage-duplex.md#output).
