# Appendix

## Concepts and terminology

| Term | Definition |
|:-----|:-----------|
| Background model | Per-site repeat-length distributions and thresholds generated from a panel of normals. |
| Microsatellite | A genomic locus containing a short, repeated sequence motif. |
| Repeat-length distribution | The frequencies of observed repeat lengths at one microsatellite site. |
| MSI-H | Run-level status assigned when the MSI score is at or above the configured aggregate threshold. |
| MSS | Run-level status assigned when the MSI score is below the configured aggregate threshold. |
| MSI score | Fraction of usable microsatellite sites classified as unstable. |
| Panel of normals | A set of normal BAMs used to estimate baseline repeat-length distributions. |
| Per-site threshold | Metric value at or above which a site is classified as unstable. |
| Usable site | Site that passes the configured read-depth and read-filter requirements. |

## Bundled resources

The container stores the following ready-to-use resources under */resources*.
The command-specific usage pages identify where to use each resource.

| Resource | Purpose |
|:---------|:--------|
| *background_thresholds_99sites_ponv1_for_TO_WGS.json* | WGS background model for [tumor-only analysis](usage-tumor-only.md). It contains repeat-length distributions from a panel of normals and per-site thresholds for 99 sites. |
| *microsatellite_sites_with_thresholds_167sites_TN_WGS.tsv* | WGS site list for [tumor-normal analysis](usage-tumor-normal.md). It contains 167 sites and a per-site `t-statistic` threshold for each site. |

## Performance and cost notes

Caller runtime scales with the number of microsatellite sites, BAM read depth, and thread count.
Background-generation runtime also scales with the number of panel samples.
Resource requirements have not been characterized across all supported data types and parameter combinations.
