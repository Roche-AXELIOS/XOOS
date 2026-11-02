# Interpreting SBX Reads

## Overview

SBX sequencing produces two types of reads depending on the adapter design:

- **Duplex reads** (SBX-D) capture both strands of the same DNA molecule in a single raw read, connected by a hairpin adapter. After demultiplexing, the two strands are aligned and collapsed into a consensus sequence.
- **Simplex reads** (SBX-S) capture a single strand of each DNA molecule. Each read represents one pass of one strand, with no intramolecular consensus.

This tutorial explains how to interpret the key components of SBX reads. It is split into two chemistry-specific guides:

- **[Interpreting SBX Duplex (SBX-D) Reads](interpreting-sbx-d-reads.md)** — duplex read structure, the YC tag, and duplex quality scores.
- **[Interpreting SBX Simplex (SBX-S) Reads](interpreting-sbx-s-reads.md)** — simplex read structure and simplex quality scores.

Both guides cover:

1. **Read structure** — how reads are organized
2. **Read names** — structured identifiers encoding run, molecule, and sample information
3. **Quality scores** — how base quality is represented at each processing stage

The **YC tag** (a duplex-specific encoding that allows lossless reconstruction of the original R1 and R2 from the consensus) is covered only in the SBX-D guide.

---

## Inter-molecular consensus reads

The consensus described above is **intramolecular** — it combines the two strands of a *single* DNA molecule (duplex) into one read. This is distinct from **inter-molecular consensus**, which combines *multiple* original molecules that map to the same location into a single, higher-accuracy read.

Inter-molecular consensus reads are produced by the [Read Collapser](../read_collapser/README.md) module after alignment. The Read Collapser clusters reads by genomic position (and, when present, UMI), then collapses each cluster into a consensus sequence. This suppresses random sequencing errors and DNA-damage artifacts, and is used for applications such as ctDNA/cfDNA profiling where ultra-high sensitivity is required.

Because these reads are generated post-alignment rather than being an inherent property of the raw SBX read structure, they are not covered by the chemistry-specific guides above. See the [Read Collapser documentation](../read_collapser/README.md) for how inter-molecular consensus reads are generated, the available presets, and the resulting output files.

---

## Read Names

Both simplex and duplex read names follow the same format:

```text
{base_prefix}:{cycle_id}:{b64_uid}:{read_bitflag}|{sid_id}
```

If UMIs are present, the format extends to:

```text
{base_prefix}:{cycle_id}:{b64_uid}:{read_bitflag}|{sid_id}:{umi5p}:{umi3p}
```

### Fields

| Field | Description | Example |
|-------|-------------|---------|
| `base_prefix` | `{date}:{sequencer}:Q{queue}:R{run}` — identifies the sequencing run | `20250115:SBX-HTP01:Q1:R3` |
| `cycle_id` | Full bright cycle ID — identifies the sequencing cycle | `42` |
| `b64_uid` | Base-64 encoded unique molecule identifier (up to 6 characters) | `0A3xKf` |
| `read_bitflag` | 8-bit field (printed as decimal 0–255) indicating SID/UMI classification and threading state | `12` |
| `sid_id` | Sample identifier index (1-based integer), or `*` if unclassified | `1` |
| `umi5p` | 5′ UMI value (integer if classified, `*` otherwise) | `127` |
| `umi3p` | 3′ UMI value (integer if classified, `*` otherwise) | `84` |

### Base Prefix

The base prefix uniquely identifies the sequencing run:

```text
{date}:{sequencer}:Q{queue_num}:R{run_num}
```

- **date** — formatted as `YYYYMMDD` in the sequencing system's time zone
- **sequencer** — instrument name
- **queue_num** — run queue number (distinguishes data from different queues)
- **run_num** — run number within the queue

### read_bitflag

The `read_bitflag` is a single byte encoded as a decimal integer, following the SBX read specification. The lower 4 bits indicate which SIDs and UMIs were classified during demultiplexing; bits 7 and 6 encode threading state:

| Bit | Mask | Meaning |
|-----|------|---------|
| 7 (MSB) | `0b10000000` (128) | Read started in a threaded state |
| 6 | `0b01000000` (64) | Read ended in a threaded state |
| 3 | `0b00001000` (8) | 5′ SID classified |
| 2 | `0b00000100` (4) | 3′ SID classified |
| 1 | `0b00000010` (2) | 5′ UMI classified |
| 0 (LSB) | `0b00000001` (1) | 3′ UMI classified |

Demultiplexed duplex reads always have both SIDs classified, so the SID/UMI bits contribute at least `12` (`0b1100`). For simplex reads, a full-length read with both SIDs contributes at least `12`, while a partial read with only the 5′ SID contributes `8` (`0b1000`).

{% hint style="info" %}
For RDB input, the threading-state bits (7 and 6) are populated from the upstream `molecular_traces_read_state` annotation. When the input predates that annotation, these bits default to 0.
{% endhint %}

For chemistry-specific read name examples, see the
[SBX-D](interpreting-sbx-d-reads.md#read-name-examples) and
[SBX-S](interpreting-sbx-s-reads.md#read-name-examples) guides.
