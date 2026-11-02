# Pan-Genome Consensus Caller

## Guide contents

- [Overview](#overview): what the module does, when to use it, example
  applications, and system requirements.
- [Usage](docs/usage.md): the main command entry point and a link to the
  command page.
  - [Consensus correction](docs/usage-pangenome-consensus-caller.md) —
    running the tool, inputs, outputs, tips, and CLI options.
- [File formats](docs/file-formats.md): the corrected output BAM.
- [Troubleshooting](docs/troubleshooting.md): common issues and resolutions.
- [Appendix](docs/appendix.md): concepts and terminology, correction logic,
  and performance notes.

## Overview

The Pan-Genome Consensus Caller corrects duplex consensus base calls using a pangenome graph alignment.
During duplex consensus, if the R1 and R2 bases disagree at a position, the R1 base is typically chosen. However, a pangenome alignment of the consensus may reveal that selecting the R2 base results in a better alignment.
The module identifies these cases and updates the consensus to use the base supported by the graph.

At a glance: **Giraffe-aligned duplex BAM** → **corrected BAM**.

To use the module, input BAM files must be aligned with [Giraffe](https://github.com/vgteam/vg) using the `--add-graph-aln` parameter, which records the graph alignment in each read's `GR` tag.
For the meaning of duplex consensus, discordant base, pangenome alignment, and the `GR` and `YC` tags, see [Concepts and terminology](docs/appendix.md#concepts-and-terminology).

{% hint style="warning" %}
**Warning:** This module has undergone limited testing.
It is intended for use only with the 7 BAM files provided with the Early Access webinar release.
{% endhint %}

### Example applications

Expand the example below for a representative scenario, its objective, the command to run, and a link to the matching usage section:

<details>

<summary>Duplex consensus correction</summary>

- **Data type:** duplex reads aligned with Giraffe using `--add-graph-aln`.
- **Objective:** correct discordant duplex consensus bases using the pangenome
  graph alignment.
- **Usage:** To run this application, go to [Consensus correction](docs/usage-pangenome-consensus-caller.md).
- **Command:** `pangenome_consensus_caller -b no_repeat.bed -t 16 < input.bam > output.bam`

</details>

### System requirements

|        | Requirements                                        |
|:--------|:-----------------------------------------------------|
| CPU    | Scales with the number of threads set by `-t`.      |
| Memory | Scales with the number of threads set by `-t`.      |
