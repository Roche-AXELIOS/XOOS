#!/usr/bin/env bash
#
# Reproduce the hg38 tumor-only WGS blocklist BED for copy_number_caller.
#
# The output is passed to `copy_number_caller tumor-only-wgs` via `--blocklist-bed`.
# It merges three sources of regions that produce spurious copy-number segments
# in tumor-only WGS:
#
#   1. gnomAD-SV v4.1 common germline CNVs (PASS DEL/DUP, global AF > 0.1%)
#   2. hg38 blocklist (the module resource resources/hg38_blocklist.bed)
#   3. the HLA/MHC region (chr6:29,500,000-33,500,000), whose extreme
#      population copy-number and mappability variability yields false
#      deep-deletion segments that bias the purity/ploidy search.
#
# Reference build: GRCh38 (no-alt), "chr"-prefixed contigs.
#
# This formalizes a merge that was previously performed in an internal R&D workflow.
# It is provided here so the blocklist resource can be reproduced from public inputs.
#
# Required tools: wget, bcftools, bedtools, sort, cut.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPT_DIR
MODULE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
readonly MODULE_DIR

# ---------------------------------------------------------------------------
# Defaults (override via flags or environment variables)
# ---------------------------------------------------------------------------
GNOMAD_URL="${GNOMAD_URL:-https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz}"
GNOMAD_VCF="${GNOMAD_VCF:-}"
BLOCKLIST_BED="${BLOCKLIST_BED:-${MODULE_DIR}/resources/hg38_blocklist.bed}"
OUTPUT_BED="${OUTPUT_BED:-${MODULE_DIR}/resources/hg38_tumor_only_wgs_blocklist.bed}"
WORK_DIR="${WORK_DIR:-}"

# HLA/MHC region to mask (hg38 coordinates, 0-based half-open).
readonly HLA_CHROM="chr6"
readonly HLA_START="29500000"
readonly HLA_END="33500000"

usage() {
    cat <<EOF
Usage: ${0##*/} [options]

Build the hg38 tumor-only WGS blocklist BED by merging gnomAD-SV common CNVs,
the hg38 blocklist, and the HLA/MHC region.

Options:
  -o, --output PATH       Output BED path
                          (default: resources/hg38_tumor_only_wgs_blocklist.bed)
  -b, --blocklist PATH    hg38 blocklist BED input
                          (default: resources/hg38_blocklist.bed; requires
                          'dvc pull resources/hg38_blocklist.bed.dvc')
  -g, --gnomad-vcf PATH   Use an already-downloaded gnomAD-SV sites VCF instead
                          of downloading it
  -u, --gnomad-url URL    gnomAD-SV sites VCF download URL
  -w, --work-dir PATH     Working directory for downloads and intermediates
                          (default: a temp dir, removed on exit)
  -h, --help              Show this help and exit

Each option has a matching environment variable: OUTPUT_BED, BLOCKLIST_BED,
GNOMAD_VCF, GNOMAD_URL, WORK_DIR.
EOF
}

log() { printf '[%s] %s\n' "${0##*/}" "$*" >&2; }
die() { printf '[%s] ERROR: %s\n' "${0##*/}" "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output)      OUTPUT_BED="$2"; shift 2 ;;
        -b|--blocklist)   BLOCKLIST_BED="$2"; shift 2 ;;
        -g|--gnomad-vcf)  GNOMAD_VCF="$2"; shift 2 ;;
        -u|--gnomad-url)  GNOMAD_URL="$2"; shift 2 ;;
        -w|--work-dir)    WORK_DIR="$2"; shift 2 ;;
        -h|--help)        usage; exit 0 ;;
        *)                usage >&2; die "unknown argument: $1" ;;
    esac
done

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
for tool in bcftools bedtools sort cut; do
    command -v "${tool}" >/dev/null 2>&1 || die "required tool not on PATH: ${tool}"
done

if [[ ! -f "${BLOCKLIST_BED}" ]]; then
    die "hg38 blocklist BED not found: ${BLOCKLIST_BED}
Run 'dvc pull resources/hg38_blocklist.bed.dvc' from the module root first, or
pass --blocklist PATH."
fi

# ---------------------------------------------------------------------------
# Working directory (auto-created temp dir is cleaned up on exit)
# ---------------------------------------------------------------------------
cleanup_work_dir=false
if [[ -z "${WORK_DIR}" ]]; then
    WORK_DIR="$(mktemp -d)"
    cleanup_work_dir=true
else
    mkdir -p "${WORK_DIR}"
fi
cleanup() { [[ "${cleanup_work_dir}" == true ]] && rm -rf "${WORK_DIR}"; }
trap cleanup EXIT

# ---------------------------------------------------------------------------
# 1. Obtain the gnomAD-SV sites VCF
# ---------------------------------------------------------------------------
if [[ -z "${GNOMAD_VCF}" ]]; then
    command -v wget >/dev/null 2>&1 || die "wget required to download gnomAD-SV (or pass --gnomad-vcf)"
    GNOMAD_VCF="${WORK_DIR}/gnomad.v4.1.sv.sites.vcf.gz"
    if [[ ! -f "${GNOMAD_VCF}" ]]; then
        log "Downloading gnomAD-SV v4.1 sites VCF (~1.7 GB) from ${GNOMAD_URL}"
        wget -q --show-progress -O "${GNOMAD_VCF}" "${GNOMAD_URL}"
    else
        log "SKIP download: ${GNOMAD_VCF} already exists"
    fi
else
    [[ -f "${GNOMAD_VCF}" ]] || die "gnomAD-SV VCF not found: ${GNOMAD_VCF}"
    log "Using provided gnomAD-SV VCF: ${GNOMAD_VCF}"
fi

# ---------------------------------------------------------------------------
# 2. Filter gnomAD-SV to common (AF > 0.1%) PASS DEL/DUP events.
#    Columns: chrom, start (0-based), end
# ---------------------------------------------------------------------------
gnomad_bed="${WORK_DIR}/gnomad_sv_v4.1_common_cnv_af0.1pct.bed"
log "Filtering gnomAD-SV to common (AF>0.1%) DEL/DUP events"
bcftools query \
        -f '%CHROM\t%POS0\t%INFO/END\n' \
        -i 'FILTER="PASS" && (INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP") && INFO/AF>0.001' \
        "${GNOMAD_VCF}" \
    | sort -k1,1V -k2,2n \
    > "${gnomad_bed}"
log "gnomAD-SV common CNV regions: $(wc -l < "${gnomad_bed}")"

# ---------------------------------------------------------------------------
# 3. Write the HLA/MHC masking BED (hg38).
# ---------------------------------------------------------------------------
hla_bed="${WORK_DIR}/hg38_hla.bed"
printf '%s\t%s\t%s\n' "${HLA_CHROM}" "${HLA_START}" "${HLA_END}" > "${hla_bed}"

# ---------------------------------------------------------------------------
# 4. Merge the three sources into a single sorted, non-overlapping BED3.
# ---------------------------------------------------------------------------
log "Merging gnomAD-SV common CNVs, hg38 blocklist, and HLA/MHC region"
mkdir -p "$(dirname "${OUTPUT_BED}")"
cat \
        <(cut -f1-3 "${gnomad_bed}") \
        <(cut -f1-3 "${BLOCKLIST_BED}") \
        <(cut -f1-3 "${hla_bed}") \
    | sort -k1,1V -k2,2n \
    | bedtools merge \
    > "${OUTPUT_BED}"

log "Wrote ${OUTPUT_BED} ($(wc -l < "${OUTPUT_BED}") merged regions)"
