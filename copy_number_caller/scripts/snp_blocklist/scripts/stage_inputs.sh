#!/usr/bin/env bash
#
# Stage the SNP blocklist inputs so config/blocklist.yaml is runnable.
#
# Downloads the two input sets needed to reproduce the reference blocklist:
#   1. The 50 raw mutect2 tumour-normal VCFs (+ .tbi), from the S3 data
#      catalog. These are the panel referenced by config/samples.tsv.
#   2. The public af-only gnomAD hg38 VCF (GATK best-practices), from GCS.
#
# The VCFs live in the S3 data catalog so the workflow no longer depends on
# any site-specific filesystem path. They are pulled to inputs/vcfs/ where
# config/samples.tsv expects them (relative inputs/vcfs/<file> paths).
#
# Prerequisites:
#   * aws CLI with credentials that can read the data-catalog bucket
#     (the bucket enforces SSE-KMS; a standard GetObject role is sufficient).
#   * curl for the public gnomAD download.
#
# Usage (from scripts/snp_blocklist):
#   scripts/stage_inputs.sh
#   snakemake -s snp_blocklist.smk --configfile config/blocklist.yaml \
#       --use-singularity --cores 8
set -euo pipefail

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# Directory of this script, then the workflow root (its parent).
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
workflow_dir=$(dirname "${script_dir}")
cd "${workflow_dir}"

sample_sheet="config/samples.tsv"
vcf_dir="inputs/vcfs"
s3_vcf_prefix="s3://prod-data-catalog-20221116223806161200000001/cnc_resources/snp_blocklist_bed/inputs/vcfs"
gnomad_dir="inputs/gnomad_af"
gnomad_vcf="${gnomad_dir}/af-only-gnomad.hg38.vcf.gz"
gnomad_url="https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"

# --- 1. Download the raw mutect2 T-N VCFs from the S3 data catalog ----------
mkdir -p "${vcf_dir}"
if ! command -v aws >/dev/null 2>&1; then
    log "ERROR: aws CLI not found; needed to download VCFs from ${s3_vcf_prefix}"
    exit 1
fi
log "Downloading raw mutect2 T-N VCFs (~144 GiB) from ${s3_vcf_prefix}"
aws s3 sync "${s3_vcf_prefix}/" "${vcf_dir}/" \
    --exclude "*" --include "*.vcf.gz" --include "*.vcf.gz.tbi"

# Verify every VCF referenced in the sample sheet is now present locally.
# Sheet VCF paths are resolved relative to the sheet's own directory (matching
# snp_blocklist.smk), so a path such as ../inputs/vcfs/<file> lands at the
# workflow-root inputs/vcfs/ this script just synced.
log "Verifying VCFs from ${sample_sheet}"
sheet_dir=$(cd "$(dirname "${sample_sheet}")" && pwd)
missing=0
# Columns: sample <TAB> vcf <TAB> normal_sample (skip header).
while IFS=$'\t' read -r case vcf _; do
    [[ -n "${case}" ]] || continue
    resolved="${vcf}"
    [[ "${resolved}" = /* ]] || resolved="${sheet_dir}/${vcf}"
    if [[ ! -e "${resolved}" ]]; then
        log "WARNING: missing VCF for ${case}: ${resolved}"
        missing=$((missing + 1))
    fi
done < <(tail -n +2 "${sample_sheet}")
if (( missing > 0 )); then
    log "WARNING: ${missing} VCF(s) missing after sync (check S3 access)"
else
    log "All VCFs present"
fi

# --- 2. Download af-only gnomAD VCF ----------------------------------------
mkdir -p "${gnomad_dir}"
if [[ -e "${gnomad_vcf}" ]]; then
    log "gnomAD AF VCF already present: ${gnomad_vcf}"
else
    log "Downloading af-only gnomAD hg38 VCF (~3 GB)"
    curl -fSL -o "${gnomad_vcf}" "${gnomad_url}"
    curl -fSL -o "${gnomad_vcf}.tbi" "${gnomad_url}.tbi"
fi
log "gnomAD AF VCF ready: ${gnomad_vcf}"

log "Done. Run: snakemake -s snp_blocklist.smk \\"
log "    --configfile config/blocklist.yaml --use-singularity --cores 8"
