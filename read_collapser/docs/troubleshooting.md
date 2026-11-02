# Troubleshooting

## High memory usage or crash on large datasets

**Symptom:** Read Collapser consumes a large amount of memory and crashes,
especially on devices with less than 16 GiB of RAM.

**Cause:** The default batch size can require significant memory on large
datasets.

**Resolution:** Use a smaller batch size or increase the available memory. See
the [Performance and cost notes](appendix.md#performance-and-cost-notes) in the
Appendix for details on batching.

## Invalid input error: unsupported format

**Symptom:** The command fails with an invalid input error.

**Cause:** Read Collapser only accepts BAM input. SAM is not accepted and raises
an error.

**Resolution:** Provide a sorted BAM file as input. Convert any SAM files to BAM
using standard bioinformatics tools.

## Invalid input error: bad index

**Symptom:** Read Collapser accepts a valid BAM file but raises a bad access
error from an invalid index.

**Cause:** The BAM index does not match the input file.

**Resolution:** Confirm the index is appropriate. Reindex the BAM input with
standard bioinformatics tools if needed.

## Calculated consensus matrix width exceeds maximum possible width

**Symptom:** Consensus fails with the error message
`Calculated consensus matrix width exceeds maximum possible width`.

**Cause:** The calculated width of the consensus matrix is larger than the
maximum allowed width. This can happen if the alignment positions are incorrect
or if there is a clustering error, and likely indicates a problem upstream in the
clustering or duplex read decoding (if enabled) process.

**Resolution:** File a bug report with the relevant details.

## Read name parsing error: does not contain '|'

**Symptom:** The command fails with `<read_name> does not contain '|'`.

**Cause:** Read Collapser parses read names to determine adapter information
(full/partial read status and UMI presence). This error means the input BAM
contains reads whose names do not follow the expected SBX read name format. This
typically indicates the BAM was not produced by the SBX demux pipeline, or the
reads were modified by an upstream tool that stripped the read name metadata.

**Resolution:** Investigate why the read names are not in the expected format. If
you are certain your data consists entirely of full-length reads with no UMIs,
you can use the hidden option `--ignore-read-name-parsing-errors` to bypass this
check. Unparseable reads are treated as full reads with both SIDs present and no
UMIs. Do not use this option if your data is expected to contain UMIs or partial
reads, because it produces incorrect results without warning.

If you have issues or concerns not documented here, please create an issue on
[GitHub](https://github.com/Roche-AXELIOS/XOOS/issues) or reach out to your local
Roche affiliate.
