# Build an improved panel-of-normals SNP blocklist (Bayes-factor model check).
#
# Origin: ported from an internal analysis prototype
# (libs/snakemake/10_improved-snp-pon.smk) into the XOOS Copy Number Caller.
# This is the self-contained blocklist path only: it starts from per-normal
# filtered SNP VCFs and produces a single snp_blocklist.bed.
#
# Rather than calling a genotype, each site is scored by the log10 Bayes factor
# comparing a weighted mixture of clean germline genotype models
# (hom_ref / het / hom_alt) against a diffuse Beta(1, 1) alternative. A strongly
# negative CLEAN_LOG10BF means the allele counts do not fit any clean germline
# genotype (e.g. mapping/allele-balance skew), so the site is a blocklist
# candidate.
#
#   Step 1 (classify_genotypes_bayes_factor): annotate each site with a
#     CLEAN_LOG10BF FORMAT field (Float).
#   Step 2 (slim_normal_genotypes): strip each scored VCF down to CLEAN_LOG10BF
#     and write as BCF (makes the merge fast).
#   Step 3 (merge_normal_genotypes): merge the slimmed per-sample BCFs into a
#     multi-sample VCF for cross-normal comparison.
#   Step 4 (extract_blocklist): write a BED of sites that are skewed
#     (CLEAN_LOG10BF < bf_threshold) in at least min_skewed_fraction of
#     >= min_informative_samples informative normals.
#
# The workflow accepts inputs at either of two boundaries, controlled by the
# prefilter flag:
#   * prefilter=False (default): inputs are per-normal filtered SNP VCFs
#     (normal_dp15) produced upstream; they are scored directly.
#   * prefilter=True: inputs are raw somatic (mutect2) tumour-normal VCFs; the
#     workflow first runs the pre-filter (sort_vcf + filter_normal_snps)
#     to derive the per-normal normal_dp15 VCFs, then scores them. This makes it
#     possible to reproduce snp_blocklist.bed from raw VCFs, and to build a fresh
#     blocklist from any new set of raw T-N VCFs.
#
# Inputs can be supplied two ways (exactly one is required):
#   1. Directory + glob (ad hoc): point input_dir at a directory of VCFs; the
#      sample name is derived from each filename.
#   2. Sample sheet (reproducible): a TSV with a `sample<TAB>vcf` header (plus an
#      optional `normal_sample` column) naming each sample and its VCF path
#      explicitly. Use this to pin a fixed, committed panel (e.g. the 50 DLS
#      normals) independent of on-disk filenames. See config/blocklist.yaml for
#      the shipped example.
#
# Usage:
#   # ad hoc directory + glob over already-filtered per-normal VCFs
#   snakemake -s snp_blocklist.smk \
#       --config input_dir=/path/to/normal_vcfs output_dir=results \
#       --cores 8
#
#   # reproducible panel from raw T-N VCFs (prefilter on) via a
#   # committed config + sample sheet
#   snakemake -s snp_blocklist.smk \
#       --configfile config/blocklist.yaml --cores 8
#
# Configuration (override with --config key=value):
#   input_dir               Directory containing input VCFs
#                                       [required unless sample_sheet]
#   input_glob              Glob to match VCFs within input_dir
#                                                 [*.normal_dp15.vcf.gz]
#   sample_sheet            TSV (sample<TAB>vcf[<TAB>normal_sample] header) of
#                                       samples [required unless input_dir]
#   output_dir              Output directory                       [outputs]
#   prefilter               Treat inputs as raw T-N VCFs and run the
#                                       normal_dp15 pre-filter first  [False]
#   population_vcf          Population AF VCF for the pre-filter
#                                       [required when prefilter=True]
#   min_normal_depth        Minimum normal-sample depth in pre-filter [15]
#   normal_sample_suffix    Normal column = sample + this suffix    [_N]
#   bf_threshold            CLEAN_LOG10BF below this is "skewed"    [-0.5]
#   min_skewed_fraction     Fraction of informative normals skewed [0.25]
#   min_informative_samples Minimum informative normals per site   [10]
#   weight_hom_ref          hom_ref mixture weight                 [0.72]
#   weight_het              het mixture weight                     [0.19]
#   weight_hom_alt          hom_alt mixture weight                 [0.09]
#   classifier              Path to classify_genotypes_bayes_factor.py
#                                       [alongside this workflow]
#   python                  Python interpreter for the classifier
#                                       (must have pysam)          [python3]

import csv
import glob
import os

# --- Configuration with defaults ---------------------------------------------
WORKFLOW_DIR = os.path.dirname(workflow.snakefile)

INPUT_DIR = config.get("input_dir")
SAMPLE_SHEET = config.get("sample_sheet")
if bool(INPUT_DIR) == bool(SAMPLE_SHEET):
    raise ValueError(
        "Provide exactly one input mode: either config['input_dir'] "
        "(a directory of per-normal filtered SNP VCFs) or "
        "config['sample_sheet'] (a TSV with a 'sample<TAB>vcf' header). "
        "Got "
        + (
            "both"
            if INPUT_DIR and SAMPLE_SHEET
            else "neither"
        )
        + "."
    )

INPUT_GLOB = config.get("input_glob", "*.normal_dp15.vcf.gz")
OUTPUT_DIR = config.get("output_dir", "outputs")
LOGS_DIR = os.path.join(OUTPUT_DIR, "logs")

# When prefilter is enabled the discovered VCFs are raw somatic (mutect2) T-N
# VCFs, and the workflow runs the pre-filter (sort + biallelic-SNP /
# normal-depth / population-overlap filter, keeping only the normal sample)
# to derive the per-normal normal_dp15 VCFs before scoring. When disabled the
# discovered VCFs are assumed to already be per-normal filtered SNP VCFs.
PREFILTER = bool(config.get("prefilter", False))

# Population allele-frequency VCF (e.g. af-only-gnomAD hg38), used both to
# restrict the pre-filter to known SNP sites and, optionally, to derive the
# genotype priors. Required when prefilter is enabled.
POPULATION_VCF = config.get("population_vcf") or ""
if PREFILTER and not POPULATION_VCF:
    raise ValueError(
        "config['prefilter'] is set but config['population_vcf'] is not. "
        "Provide the population allele-frequency VCF (e.g. "
        "af-only-gnomad.hg38.vcf.gz) used to restrict the SNP pre-filter."
    )

# Minimum normal-sample depth required to keep a site in the pre-filter.
MIN_NORMAL_DEPTH = config.get("min_normal_depth", 15)

# Suffix appended to a sample name to form its normal-sample column name in the
# raw T-N VCF (e.g. sample DLS_7 -> normal column DLS_7_N). Used only in
# prefilter mode, and only when a sample sheet does not give an explicit
# normal_sample column.
NORMAL_SAMPLE_SUFFIX = config.get("normal_sample_suffix", "_N")

BF_THRESHOLD = config.get("bf_threshold", -0.5)
MIN_SKEWED_FRACTION = config.get("min_skewed_fraction", 0.25)
MIN_INFORMATIVE_SAMPLES = config.get("min_informative_samples", 10)

WEIGHT_HOM_REF = config.get("weight_hom_ref", 0.72)
WEIGHT_HET = config.get("weight_het", 0.19)
WEIGHT_HOM_ALT = config.get("weight_hom_alt", 0.09)

CLASSIFIER = config.get(
    "classifier",
    os.path.join(WORKFLOW_DIR, "classify_genotypes_bayes_factor.py"),
)

# Python interpreter used to run the classifier (must have pysam available).
PYTHON = config.get("python", "python3")

BCFTOOLS_CONTAINER = "docker://quay.io/biocontainers/bcftools:1.20--h8b25389_0"


# --- Derive the sample set --------------------------------------------------
def _discover_from_glob() -> dict[str, str]:
    """Map sample name -> input VCF path for every VCF matching the glob.

    The sample name is the input filename with the standard VCF suffixes and
    the ``.normal_dp15`` tag stripped, so ``DLS_7.normal_dp15.vcf.gz`` yields
    the sample ``DLS_7``.
    """
    matches = sorted(glob.glob(os.path.join(INPUT_DIR, INPUT_GLOB)))
    if not matches:
        raise ValueError(
            f"No input VCFs matched {INPUT_GLOB!r} in {INPUT_DIR!r}. "
            "Check config['input_dir'] and config['input_glob']."
        )
    samples: dict[str, str] = {}
    for path in matches:
        name = os.path.basename(path)
        for suffix in (".vcf.gz", ".vcf.bgz", ".vcf"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        if name.endswith(".normal_dp15"):
            name = name[: -len(".normal_dp15")]
        if name in samples:
            raise ValueError(
                f"Duplicate sample name {name!r} derived from "
                f"{path!r} and {samples[name]!r}."
            )
        samples[name] = path
    return samples


def _discover_from_sample_sheet() -> tuple[dict[str, str], dict[str, str]]:
    """Read a TSV sample sheet into (sample -> VCF, sample -> normal column).

    The sheet is tab-separated with a ``sample<TAB>vcf`` header and an optional
    third ``normal_sample`` column. Blank lines and ``#`` comment lines are
    ignored, and ``vcf`` paths are resolved relative to the sample sheet's own
    directory when not absolute. This makes the panel reproducible from a
    committed file regardless of the working directory.

    The ``normal_sample`` column names the normal-sample VCF column and is only
    used in prefilter mode; when omitted it defaults to
    ``sample + normal_sample_suffix``.
    """
    sheet_dir = os.path.dirname(os.path.abspath(SAMPLE_SHEET))
    allowed = (["sample", "vcf"], ["sample", "vcf", "normal_sample"])
    samples: dict[str, str] = {}
    normals: dict[str, str] = {}
    with open(SAMPLE_SHEET, newline="") as handle:
        rows = [
            line
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        ]
        reader = csv.DictReader(rows, delimiter="\t")
        if reader.fieldnames not in allowed:
            raise ValueError(
                f"Sample sheet {SAMPLE_SHEET!r} must have a "
                f"'sample<TAB>vcf' header (with an optional "
                f"'normal_sample' column), got {reader.fieldnames!r}."
            )
        for row in reader:
            name = (row["sample"] or "").strip()
            vcf = (row["vcf"] or "").strip()
            if not name or not vcf:
                raise ValueError(
                    f"Sample sheet {SAMPLE_SHEET!r} has a row with an "
                    f"empty sample or vcf field: {row!r}."
                )
            if not os.path.isabs(vcf):
                vcf = os.path.normpath(os.path.join(sheet_dir, vcf))
            if name in samples:
                raise ValueError(
                    f"Duplicate sample name {name!r} in sample sheet "
                    f"{SAMPLE_SHEET!r}."
                )
            samples[name] = vcf
            normal = (row.get("normal_sample") or "").strip()
            if normal:
                normals[name] = normal
    if not samples:
        raise ValueError(
            f"Sample sheet {SAMPLE_SHEET!r} lists no samples."
        )
    return samples, normals


if SAMPLE_SHEET:
    SAMPLE_VCFS, SAMPLE_NORMALS = _discover_from_sample_sheet()
else:
    SAMPLE_VCFS, SAMPLE_NORMALS = _discover_from_glob(), {}
CASES = sorted(SAMPLE_VCFS)


def _normal_sample_name(case: str) -> str:
    """Return the normal-sample column name for a case in prefilter mode.

    Uses the explicit ``normal_sample`` from the sample sheet when given,
    otherwise falls back to ``case + normal_sample_suffix``.
    """
    return SAMPLE_NORMALS.get(case, case + NORMAL_SAMPLE_SUFFIX)


rule all:
    input:
        f"{OUTPUT_DIR}/extract_blocklist/snp_blocklist.bed",


def _scorer_input(wildcards) -> str:
    """Per-normal filtered SNP VCF fed to the Bayes-factor scorer.

    In prefilter mode this is the filter_normal_snps output derived from the
    raw T-N VCF; otherwise it is the discovered per-normal VCF used directly.
    """
    if PREFILTER:
        return (
            f"{OUTPUT_DIR}/filter_normal_snps/"
            f"{wildcards.case}.normal_dp15.vcf.gz"
        )
    return SAMPLE_VCFS[wildcards.case]


# --- Pre-filter: raw somatic T-N VCF -> per-normal filtered SNP VCF ----------
# These two rules run only in prefilter mode. They reproduce the
# pre-filter that turns a raw mutect2 tumour-normal VCF into the per-normal
# normal_dp15 VCF the scorer expects.
rule sort_vcf:
    """Sort and index a raw somatic (mutect2) T-N VCF with bcftools."""
    input:
        vcf=lambda wildcards: SAMPLE_VCFS[wildcards.case],
    output:
        vcf=f"{OUTPUT_DIR}/sort_vcf/{{case}}_mutect2.vcf.gz",
        tbi=f"{OUTPUT_DIR}/sort_vcf/{{case}}_mutect2.vcf.gz.tbi",
    threads: 1
    log:
        stdout=f"{LOGS_DIR}/sort_vcf/{{case}}.stdout",
        stderr=f"{LOGS_DIR}/sort_vcf/{{case}}.stderr",
    container:
        BCFTOOLS_CONTAINER
    shell:
        """
        bcftools sort \\
                --output-type z \\
                --output {output.vcf} \\
                {input.vcf} \\
            2> {log.stderr} \\
            > {log.stdout}

        bcftools index --tbi {output.vcf} 2>> {log.stderr}
        """


rule filter_normal_snps:
    """Filter a sorted T-N VCF to per-normal biallelic SNPs at good depth.

    Keeps biallelic SNPs where the normal-sample depth exceeds
    min_normal_depth and the site overlaps the population VCF, retaining only
    the normal sample. This is the normal_dp15 filter.
    """
    input:
        vcf=f"{OUTPUT_DIR}/sort_vcf/{{case}}_mutect2.vcf.gz",
        tbi=f"{OUTPUT_DIR}/sort_vcf/{{case}}_mutect2.vcf.gz.tbi",
        population_vcf=POPULATION_VCF,
    output:
        vcf=f"{OUTPUT_DIR}/filter_normal_snps/{{case}}.normal_dp15.vcf.gz",
        tbi=f"{OUTPUT_DIR}/filter_normal_snps/{{case}}.normal_dp15.vcf.gz.tbi",
    params:
        normal_sample_name=lambda wildcards: _normal_sample_name(
            wildcards.case
        ),
        min_normal_depth=MIN_NORMAL_DEPTH,
    threads: 4
    log:
        stdout=f"{LOGS_DIR}/filter_normal_snps/{{case}}.stdout",
        stderr=f"{LOGS_DIR}/filter_normal_snps/{{case}}.stderr",
    container:
        BCFTOOLS_CONTAINER
    shell:
        """
        normal_sample_idx=$(bcftools query -l {input.vcf} \\
            | awk '$0=="{params.normal_sample_name}"{{print NR-1; exit}}')
        if [[ -z "$normal_sample_idx" ]]; then
            echo "Normal sample {params.normal_sample_name} not found in \\
{input.vcf}" > {log.stderr}
            exit 1
        fi
        bcftools view \\
                --min-alleles 2 \\
                --max-alleles 2 \\
                --types snps \\
                --targets-file {input.population_vcf} \\
                --include "FMT/DP[$normal_sample_idx] > \\
{params.min_normal_depth}" \\
                --samples {params.normal_sample_name} \\
                --output-type z \\
                --output {output.vcf} \\
                --threads {threads} \\
                {input.vcf} \\
            2> {log.stderr} \\
            > {log.stdout}

        bcftools index --tbi {output.vcf} 2>> {log.stderr}
        """


rule classify_genotypes_bayes_factor:
    """Annotate each site with the clean-genotype log10 Bayes factor.

    Scores a per-sample filtered normal VCF by the log10 Bayes factor of a
    clean-genotype Beta-Binomial mixture vs a diffuse Beta(1, 1) alternative,
    written as a new FORMAT field CLEAN_LOG10BF; strongly negative values flag
    sites that do not fit a clean germline genotype.
    """
    input:
        vcf=_scorer_input,
        script=CLASSIFIER,
    output:
        vcf=f"{OUTPUT_DIR}/genotype_classification/{{case}}.normal_dp15.genotypes.vcf.gz",
        tbi=f"{OUTPUT_DIR}/genotype_classification/{{case}}.normal_dp15.genotypes.vcf.gz.tbi",
    params:
        weight_hom_ref=WEIGHT_HOM_REF,
        weight_het=WEIGHT_HET,
        weight_hom_alt=WEIGHT_HOM_ALT,
    threads: 1
    log:
        stdout=f"{LOGS_DIR}/classify_genotypes_bayes_factor/{{case}}.stdout",
        stderr=f"{LOGS_DIR}/classify_genotypes_bayes_factor/{{case}}.stderr",
    shell:
        """
        {PYTHON} {input.script} \
                --input {input.vcf} \
                --output {output.vcf} \
                --weight-hom-ref {params.weight_hom_ref} \
                --weight-het {params.weight_het} \
                --weight-hom-alt {params.weight_hom_alt} \
            > {log.stdout} \
            2> {log.stderr}

        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """


rule slim_normal_genotypes:
    """Strip each scored VCF down to just the CLEAN_LOG10BF FORMAT field.

    Dropping INFO and all other FORMAT fields (GT/AD/...) and writing BCF makes
    the downstream bcftools merge far faster: merge then only has to decompress
    and harmonize the single field it actually needs.
    """
    input:
        vcf=f"{OUTPUT_DIR}/genotype_classification/{{case}}.normal_dp15.genotypes.vcf.gz",
        tbi=f"{OUTPUT_DIR}/genotype_classification/{{case}}.normal_dp15.genotypes.vcf.gz.tbi",
    output:
        bcf=f"{OUTPUT_DIR}/slim_normal_genotypes/{{case}}.clean_log10bf.bcf",
        csi=f"{OUTPUT_DIR}/slim_normal_genotypes/{{case}}.clean_log10bf.bcf.csi",
    threads: 1
    log:
        stdout=f"{LOGS_DIR}/slim_normal_genotypes/{{case}}.stdout",
        stderr=f"{LOGS_DIR}/slim_normal_genotypes/{{case}}.stderr",
    container:
        BCFTOOLS_CONTAINER
    shell:
        """
        bcftools annotate \
                --remove 'INFO,^FORMAT/CLEAN_LOG10BF' \
                --output-type b \
                --output {output.bcf} \
                {input.vcf} \
            2> {log.stderr} \
            > {log.stdout}

        bcftools index {output.bcf} 2>> {log.stderr}
        """


rule merge_normal_genotypes:
    """Merge the slimmed per-sample CLEAN_LOG10BF BCFs into a multi-sample VCF.

    Merging the slimmed BCFs (CLEAN_LOG10BF only) is much faster than merging
    the full genotype VCFs, and yields a matrix of CLEAN_LOG10BF values that can
    be compared across normals.
    """
    input:
        bcfs=expand(
            f"{OUTPUT_DIR}/slim_normal_genotypes/{{case}}.clean_log10bf.bcf",
            case=CASES,
        ),
        csis=expand(
            f"{OUTPUT_DIR}/slim_normal_genotypes/{{case}}.clean_log10bf.bcf.csi",
            case=CASES,
        ),
    output:
        vcf=f"{OUTPUT_DIR}/merge_normal_genotypes/all_samples.normal_dp15.genotypes.vcf.gz",
        tbi=f"{OUTPUT_DIR}/merge_normal_genotypes/all_samples.normal_dp15.genotypes.vcf.gz.tbi",
    threads: 8
    log:
        stdout=f"{LOGS_DIR}/merge_normal_genotypes/all_samples.stdout",
        stderr=f"{LOGS_DIR}/merge_normal_genotypes/all_samples.stderr",
    container:
        BCFTOOLS_CONTAINER
    shell:
        """
        bcftools merge \
                --merge none \
                --output-type z \
                --output {output.vcf} \
                --threads {threads} \
                {input.bcfs} \
            2> {log.stderr} \
            > {log.stdout}

        bcftools index --tbi {output.vcf} 2>> {log.stderr}
        """


rule extract_blocklist:
    """Write a BED of SNP sites that do not fit a clean germline genotype.

    A site is blocklisted when at least min_informative_samples normals have a
    non-missing CLEAN_LOG10BF and at least min_skewed_fraction of those
    informative normals are skewed (CLEAN_LOG10BF < bf_threshold).  Samples with
    a missing CLEAN_LOG10BF (no AD at the site) do not count as informative.
    """
    input:
        vcf=f"{OUTPUT_DIR}/merge_normal_genotypes/all_samples.normal_dp15.genotypes.vcf.gz",
        tbi=f"{OUTPUT_DIR}/merge_normal_genotypes/all_samples.normal_dp15.genotypes.vcf.gz.tbi",
    output:
        bed=f"{OUTPUT_DIR}/extract_blocklist/snp_blocklist.bed",
    params:
        bf_threshold=BF_THRESHOLD,
        min_skewed_fraction=MIN_SKEWED_FRACTION,
        min_informative_samples=MIN_INFORMATIVE_SAMPLES,
    threads: 1
    log:
        stdout=f"{LOGS_DIR}/extract_blocklist/all_samples.stdout",
        stderr=f"{LOGS_DIR}/extract_blocklist/all_samples.stderr",
    container:
        BCFTOOLS_CONTAINER
    shell:
        """
        bcftools query \
                --format '%CHROM\\t%POS0\\t%END[\\t%CLEAN_LOG10BF]\\n' \
                {input.vcf} \
            2> {log.stderr} \
            | awk \
                -v thr={params.bf_threshold} \
                -v frac={params.min_skewed_fraction} \
                -v min_info={params.min_informative_samples} '
                {{
                    informative = 0
                    skewed = 0
                    for (i = 4; i <= NF; i++) {{
                        if ($i != ".") {{
                            informative++
                            if ($i + 0 < thr) skewed++
                        }}
                    }}
                    if (informative >= min_info && skewed >= frac * informative)
                        print $1 "\\t" $2 "\\t" $3
                }}' \
            > {output.bed} \
            2>> {log.stderr}
        echo "Wrote blocklist to {output.bed}" > {log.stdout}
        """
