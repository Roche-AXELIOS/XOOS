# Read Collapser

## Guide contents

- [Overview](#overview): what Read Collapser does, example applications, and system requirements.
- [Usage](docs/usage.md): running each command.
  - [Mark duplicates](docs/usage-markdup.md): running `markdup`, inputs, outputs, tips, and full CLI options.
  - [Consensus](docs/usage-consensus.md): running `consensus`, inputs, outputs, tips, and full CLI options.
  - [Rescue secondary alignments](docs/usage-rescue-secondary.md): running `read_collapser_rescue_secondary` on pangenome-aligned data, inputs, outputs, tips, and full CLI options.
- [File formats](docs/file-formats.md): custom BAM tags, consensus FASTQ record format, and consensus auxiliary tags.
- [Metrics](docs/metrics.md): metric file schemas and definitions.
- [Troubleshooting](docs/troubleshooting.md): common issues and resolutions.
- [Appendix](docs/appendix.md): concepts, terminology, and glossary.
- [Release notes](CHANGELOG.md): changes by release.

## Overview

Read Collapser deduplicates position-sorted, aligned sequencing reads by either
marking duplicates (`markdup`) or collapsing them into error-corrected consensus
sequences (`consensus`). Both operations first cluster similar reads, then act on
each cluster to reduce noise and improve the accuracy of downstream analysis.

At a glance: **position-sorted BAM** → **duplicate-marked BAM** (`markdup`) or
**consensus FASTQ** (`consensus`).

![High-level diagram](docs/assets/high-level-overview.svg)

### Sources of duplicates

The type of duplicates present in the data depends on the library preparation,
for example:

- **SBX-D**: Duplicates arise from linear amplification and sequencing.
- **SBX-Fast**: All duplicates come from sequencing.
- **PCR-based protocols**: In addition to the above, PCR amplification introduces
  PCR duplicates.

Expand the sections below for guidance on when to enable strand-aware and
UMI-based clustering for your library preparation.

<details>

<summary>When to cluster by strand</summary>

Read strand is determined by alignment to the reference genome. Strand-aware
clustering (via the `--cluster-by-strand` flag) groups reads only if they align
to the same strand. This is required for standard SBX-D and SBX-Fast
(parent-parent duplex) protocols, where library preparation and amplification
yield Xpandomers of the same orientation from a single original molecule.
Enforcing this constraint prevents the clustering of reads from opposite strands,
resulting in more accurate duplicate marking and consensus generation.

![Strand-aware clustering groups same-orientation reads](docs/assets/cluster-by-strand-same-orientation-text-outlined.svg)

Conversely, if the library preparation and amplification processes can produce
Xpandomers of both orientations from a single original molecule, strand-aware
clustering must be omitted. Enabling it would prevent these complementary reads
from being grouped together. Because standard protocols for SBX-S (used in target
enrichment and simplex WGS) and parent-daughter duplex sequencing produce
Xpandomers in both orientations, the `--cluster-by-strand` flag must not be
provided.

![Without strand-aware clustering, both orientations are grouped](docs/assets/no-cluster-by-strand-diff-orientation-text-outlined.svg)

The built-in presets are designed with these principles in mind, but you must
also consider the specific properties of your library preparation and
amplification process when deciding whether to enable strand-aware clustering.

</details>

<details>

<summary>When to cluster by UMI</summary>

Certain adapter types allow for the tagging of reads with UMI (unique molecular
index) tags. UMIs are short sequences that are added to both ends of the original
DNA molecule during library preparation before any amplification steps. UMIs are
trimmed from the read sequence during demultiplexing and are instead encoded in
the read name. Always enable UMI clustering (`--cluster-by-umi`) when your
library preparation includes UMIs. UMI clustering uses this information to group
reads from the same original parent molecule more accurately.

Read Collapser matches UMI sequences **exactly**; it does not perform UMI
sequence correction. Correct UMI sequences upstream before
running Read Collapser (for example, during demultiplexing). UMI
matching is combined with positional clustering, so reads must also satisfy the
applicable positional criteria.

</details>

### Example applications

Expand a section below for a representative application, including its data type,
objective, clustering technique, and the command to run. Each links to the
relevant Usage page section.

<details>
<summary>Germline WGS (duplex)</summary>

- **Data type:** SBX-D, SBX-Fast
- **Objective:** To prevent artificial inflation of read depth and ensure
  accurate quantification of inherited variants by marking system duplicates and
  any residual amplification duplicates (if applicable).
- **Clustering technique:** Position + strand
- **Command:** `markdup --preset wgs-duplex` (see [Mark duplicates: Example command](docs/usage-markdup.md#example-command))

</details>

<details>
<summary>Germline WGS (simplex)</summary>

- **Data type:** SBX-S YS P-D
- **Objective:** To mark system duplicates in simplex sequencing data while
  preserving full-length read prioritization and handling partial reads.
- **Clustering technique:** Position
- **Command:** `markdup --preset wgs-simplex` (see [Mark duplicates: Example command](docs/usage-markdup.md#example-command))

</details>

<details>
<summary>Somatic WGS TN (duplex)</summary>

- **Data type:** SBX-D FFPE
- **Objective:** To accurately quantify tumor variant allele frequencies and
  avoid biased read depth by marking duplicates from bulk sequencing.
- **Clustering technique:** Position + strand
- **Command:** `markdup --preset wgs-duplex` (see [Mark duplicates: Example command](docs/usage-markdup.md#example-command))

</details>

<details>
<summary>cfDNA WGS (duplex)</summary>

- **Data type:** SBX-D ctDNA/cfDNA
- **Objective:** To achieve ultra-high sensitivity for cfDNA profiling and
  low-frequency ctDNA detection by correcting sequencing and linear amplification
  errors at the single-molecule level.
- **Clustering technique:** Position + strand
- **Command:** `consensus --preset wgs-duplex-cfdna` (see [Consensus: Example command](docs/usage-consensus.md#example-command))

</details>

<details>
<summary>Somatic TE (simplex)</summary>

- **Data type:** SBX YSU TE ctDNA
- **Objective:** To increase the accuracy and sensitivity of somatic variant
  detection by correcting PCR and sequencing errors from reads of the same
  original molecule.
- **Clustering technique:** Position + UMI
- **Command:** `consensus --preset te-simplex` (see [Consensus: Example command](docs/usage-consensus.md#example-command))

</details>

{% hint style="warning" %}
For Small Variant Caller (SVC) workflows that expect BAM input, use the
mark-duplicate BAM from `markdup`. The ctDNA/cfDNA WGS duplex preset
(`wgs-duplex-cfdna`) outputs consensus FASTQs (and optionally a cluster BAM); any
downstream "inter-consensus BAMs" are produced only after re-aligning those
consensus FASTQs, and are intended for ctDNA/cfDNA applications such as tumor
fraction estimation and cfDNA profiling, not standard SVC WGS input.
{% endhint %}

### System requirements

|           | Requirement                                   |
|:----------|:---------------------------------------------|
| CPU       | Modern Processor with at least 6 cores.      |
| Memory    | At least 16 GiB.                             |
| Software  | Docker or equivalent containerization environment must be installed. |
