# File formats

This page describes the tool-specific formats that MSI Caller reads or writes.
Command instructions are available on the [`tumor-only`](usage-tumor-only.md), [`tumor-normal`](usage-tumor-normal.md), and [`background_generator`](usage-background-generator.md) pages.
The run-level quality-control metric schema is on the [Metrics](metrics.md) page.

## Microsatellite sites TSV

Applies to: `tumor-normal`, `background_generator`.

The sites file is a tab-separated file with a header row.
The six columns below are required, but their order does not matter and additional columns are ignored.
The table shows the canonical order.

| Field | Type | Description |
|:------|:-----|:------------|
| chromosome | String | Contig name. |
| location | Integer | Zero-based genomic start coordinate. |
| repeat_times | Integer | Reference repeat count. |
| repeat_unit_bases | String | Repeat motif. |
| left_flank_bases | String | Left flanking sequence. |
| right_flank_bases | String | Right flanking sequence. |

The `--microsatellite-sites-thresholds-tsv` form requires one additional numeric column.
Its name must match the selected metric, such as `t-statistic` or `js-distance`.

<!-- markdownlint-disable MD010 -->

```tsv
chromosome	location	repeat_times	repeat_unit_bases	left_flank_bases	right_flank_bases	t-statistic
chr2	47408384	13	T	TGTAC	AAGGA	3.628334
```

<!-- markdownlint-enable MD010 -->

Back to [`tumor-normal` input](usage-tumor-normal.md#input) or [`background_generator` input](usage-background-generator.md#input).

## Background model JSON

Applies to: `tumor-only`, `background_generator`.

The background model is a JSON array with one object per microsatellite site.

| Field | Type | Description |
|:------|:-----|:------------|
| chromosome | String | Contig name. |
| location | Nonnegative integer | Zero-based genomic start coordinate. |
| left_flank | String | Left flanking sequence. |
| reference_repeat_times | Positive integer | Reference repeat count. |
| repeat_unit | String | Repeat motif. |
| right_flank | String | Right flanking sequence. |
| start_length | Positive integer | Repeat length represented by the first distribution value. |
| distribution | Number array | Normalized background frequency for consecutive repeat lengths. |
| threshold | Object array | Metric name and threshold value pairs. |

Each `threshold` object contains the following fields:

| Field | Type | Description |
|:------|:-----|:------------|
| metric | String | Comparison metric, such as `t-statistic` or `js-distance`. |
| value | Number | Instability threshold for the metric. |

The following example contains one background site:

```json
[
  {
    "chromosome": "chr2",
    "location": 47408384,
    "left_flank": "TGTAC",
    "reference_repeat_times": 13,
    "repeat_unit": "T",
    "right_flank": "AAGGA",
    "start_length": 11,
    "distribution": [0.1, 0.7, 0.2],
    "threshold": [
      {
        "metric": "t-statistic",
        "value": 2.0
      }
    ]
  }
]
```

Back to [`tumor-only` input](usage-tumor-only.md#input) or [`background_generator` output](usage-background-generator.md#output).

## MSI JSON output

Applies to: `tumor-only`, `tumor-normal`.

The caller produces *msi_results.json* with an `inputs` object that records the run configuration and a `metrics` object that records aggregate MSI results.
The top-level `CommandLine` and `Version` fields record the command and software version that produced the report.

### Provenance fields

| Field | Type | Description |
|:------|:-----|:------------|
| CommandLine | String | Full command line used to produce the report. |
| Version | String | MSI caller version used to produce the report. |

### Input fields

| Field | Type | Description |
|:------|:-----|:------------|
| sample_name | String | Sample name written to the outputs. |
| mode | String | `TUMOR_ONLY` or `TUMOR_NORMAL`. |
| min_sites | Integer | Minimum usable-site requirement. |
| min_depth | Integer | Minimum read depth for a site. A site must meet or exceed this value. |
| mapq_cutoff | Integer | Minimum mapping quality. |
| baseq_cutoff | Integer | Minimum fallback base-call quality. |
| min_base_type | String | Minimum YC-decoded base category: `discordant`, `simplex`, or `concordant`. |
| min_query_microsatellite_length | Integer | Exclusive lower bound for query microsatellite length. |
| max_query_microsatellite_length | Integer | Exclusive upper bound for query microsatellite length. |
| trim_soft_clips | Boolean | Whether the caller trimmed soft-clipped sequence. |
| msi_classification_threshold | Number or string | MSI score threshold, or `THRESHOLD_UNSPECIFIED` when no threshold was supplied. |

In `tumor-normal` mode, the `inputs` object also contains `normal_bam`, which records the matched normal BAM path.

### Result fields

| Field | Type | Description |
|:------|:-----|:------------|
| msi_status | String | `MSI-H`, `MSS`, `THRESHOLD_UNSPECIFIED`, or `INSUFFICIENT_COVERAGE`. |
| msi_score | Number or string | Fraction of usable sites classified as unstable, or `INSUFFICIENT_COVERAGE` when no score is available. |
| total_usable_sites | Integer | Number of sites used to calculate the MSI score. |
| total_unstable_sites | Integer | Number of usable sites classified as unstable. |
| total_microsatellite_sites | Integer | Number of evaluated sites. |

The following example shows tumor-normal output:

```json
{
  "inputs": {
    "sample_name": "SAMPLE",
    "mode": "TUMOR_NORMAL",
    "normal_bam": "normal.bam",
    "min_sites": 10,
    "min_depth": 20,
    "mapq_cutoff": 0,
    "baseq_cutoff": 0,
    "min_base_type": "discordant",
    "min_query_microsatellite_length": 0,
    "max_query_microsatellite_length": 40,
    "trim_soft_clips": false,
    "msi_classification_threshold": 0.3
  },
  "metrics": {
    "msi_status": "MSI-H",
    "msi_score": 0.4,
    "total_usable_sites": 100,
    "total_unstable_sites": 40,
    "total_microsatellite_sites": 120
  },
  "CommandLine": "msi-caller --log-level info tumor-normal --tumor-bam tumor.bam --normal-bam normal.bam",
  "Version": "1.2.0"
}
```

Back to [`tumor-only` output](usage-tumor-only.md#output) or [`tumor-normal` output](usage-tumor-normal.md#output).

## MSI VCF

Applies to: `tumor-only`, `tumor-normal`.

The caller produces a VCF containing one record per configured microsatellite site.
The header includes a `##RocheCommandLine=<ID=...,Version=...,CommandLine=...>` provenance record.
The standard VCF columns describe the locus and call state.
`ID` contains the site name, `REF` contains the reference repeat sequence, and `ALT` is the symbolic allele `<MSI>`.
`FILTER` is `PASS` for a usable site or `INSUFFICIENT_COVERAGE` when read depth at that microsatellite site is below the minimum required `min_depth` parameter.   Both tumor-normal and tumor-only BAM files must meet the minimum depth requirement to `PASS`.

### Run metadata

| Field | Description |
|:------|:------------|
| MSI_STATUS | Aggregate status: `MSI-H`, `MSS`, `THRESHOLD_UNSPECIFIED`, or `INSUFFICIENT_COVERAGE`. |
| MSI_SCORE | Fraction of usable sites classified as unstable, or `INSUFFICIENT_COVERAGE` when no score is available. |
| MSI_SITES | Number of unstable sites. |
| TOTAL_USABLE_SITES | Number of sites with sufficient data for a call. |
| TOTAL_SITES | Number of evaluated sites. |
| ALT | Definition of the symbolic `<MSI>` allele. |
| THRESHOLD | Aggregate MSI classification threshold or `THRESHOLD_UNSPECIFIED`. |
| SAMPLE | Output sample name. |
| MIN_SITES | Minimum number of usable sites required for an aggregate call. |
| METRIC | Distribution-comparison metric. |
| MIN_DEPTH | Minimum read depth for a site. A site must meet or exceed this value. |
| MIN_QUERY | Exclusive lower bound for query microsatellite length. |
| MAX_QUERY | Exclusive upper bound for query microsatellite length. |
| MAPQ_CUTOFF | Minimum mapping quality. |
| BASEQ_CUTOFF | Minimum fallback base-call quality. |
| MIN_BASE_TYPE | Minimum YC-decoded base category. |
| TRIM_SOFT_CLIP | Whether soft-clipped sequence was trimmed. |
| MODE | `TUMOR_ONLY` or `TUMOR_NORMAL`. |
| FILTER | Definition of the `INSUFFICIENT_COVERAGE` record filter. |
| NORMAL_BAM | Normal BAM path in tumor-normal mode. |

### INFO fields

| Field | Description |
|:------|:------------|
| KMER | Repeat motif. |
| RPLEN | Reference repeat count. |
| THRESHOLD | Site threshold from a background model or threshold-bearing sites file. Absent with a plain sites file. |
| BGD | Background repeat-length distribution in tumor-only mode. |
| MINLEN | Minimum repeat length represented by `BGD` in tumor-only mode. |
| MAXLEN | Maximum repeat length represented by `BGD` in tumor-only mode. |

### FORMAT fields

| Field | Description |
|:------|:------------|
| SD | Tumor repeat-length distribution. |
| MINLEN | Minimum repeat length represented by `SD`. |
| MAXLEN | Maximum repeat length represented by `SD`. |
| DP | Tumor read depth used to build the distribution. |
| M_SAMPLE | Comparison metric value for the tumor site. |
| MSP | Site prediction: `U` (unstable), `S` (stable), or `N` (not predicted). |
| ND | Normal read depth in tumor-normal mode. |
| NSD | Normal repeat-length distribution in tumor-normal mode. |
| NMINLEN | Minimum repeat length represented by `NSD`. |
| NMAXLEN | Maximum repeat length represented by `NSD`. |

Back to [`tumor-only` output](usage-tumor-only.md#output) or [`tumor-normal` output](usage-tumor-normal.md#output).
