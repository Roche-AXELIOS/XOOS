# Troubleshooting

## No reads overlapping the target variants

**Symptom:** The run reports no reads overlapping the target variants.

**Cause:** All variant sites were filtered, or the BAM has no coverage in the VCF regions.

**Resolution:** Verify the BAM and VCF match the same reference genome and coordinate system.

## No reads overlapping the noise probe sites

**Symptom:** The run reports no reads overlapping the noise probe sites.

**Cause:** Noise probe generation succeeded but the BAM has no coverage at the probe positions.
This can occur with very sparse BAM files or restrictive region filters.

**Resolution:** Confirm the BAM has coverage near the query variants, and relax region filters if they are too restrictive.

## Bundled VCF not found (contamination)

**Symptom:** The `contamination` run fails because the bundled population SNP VCF cannot be found.

**Cause:** The default population SNP VCF path assumes the container image layout.

**Resolution:** When running outside Docker, provide `--query-vcf` explicitly.
