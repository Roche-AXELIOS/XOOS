#pragma once

#include <string>
#include <string_view>

#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "core/sample-type.h"

namespace xoos::svc {

// Column headers for the QC metrics TSV.
inline constexpr std::string_view kMetricColumnName = "metric_name";
inline constexpr std::string_view kMetricColumnValue = "value";
inline constexpr std::string_view kMetricColumnPercentage = "percentage";
inline constexpr std::string_view kMetricColumnDenominator = "denominator";

// Metric names (the `metric_name` column values) emitted by WriteVcfMetricsFile.
inline constexpr std::string_view kMetricQcStatus = "qc_status";
inline constexpr std::string_view kMetricPredictedSex = "predicted_sex";
inline constexpr std::string_view kMetricSampleType = "sample_type";
inline constexpr std::string_view kMetricNumInputVariants = "num_input_variants";
inline constexpr std::string_view kMetricNumOutputVariants = "num_output_variants";
inline constexpr std::string_view kMetricNumPassVariants = "num_pass_variants";
inline constexpr std::string_view kMetricNumNonPassVariants = "num_non_pass_variants";
inline constexpr std::string_view kMetricNumPassSnvs = "num_pass_snvs";
inline constexpr std::string_view kMetricNumPassInsertions = "num_pass_insertions";
inline constexpr std::string_view kMetricNumPassDeletions = "num_pass_deletions";
inline constexpr std::string_view kMetricNumPassBiallelic = "num_pass_biallelic";
inline constexpr std::string_view kMetricNumPassMultiallelic = "num_pass_multiallelic";
inline constexpr std::string_view kMetricNumPassHomozygous = "num_pass_homozygous";
inline constexpr std::string_view kMetricNumPassHeterozygous = "num_pass_heterozygous";
inline constexpr std::string_view kMetricHetHomRatio = "het_hom_ratio";
inline constexpr std::string_view kMetricNumTransitions = "num_transitions";
inline constexpr std::string_view kMetricNumTransversions = "num_transversions";
inline constexpr std::string_view kMetricTiTvRatio = "ti_tv_ratio";
// Per-FILTER metric names are formed as `<prefix><FILTER tag>`, e.g. `num_filter_ML_Score`.
inline constexpr std::string_view kMetricFilterPrefix = "num_filter_";

// Values of the `qc_status` metric.
inline constexpr std::string_view kQcStatusPass = "pass";
inline constexpr std::string_view kQcStatusFail = "fail";

/**
 * @brief Aggregated VCF-based QC metrics for the `filter_variants` output.
 *
 * Counts are expressed in terms of ACGT-only variants (ALT alleles for which both REF and ALT contain only ACGT), not
 * VCF records: a multi-allelic record with two ACGT-only ALT alleles contributes two to the variant totals and genotype
 * breakdown.
 */
struct VcfMetrics {
  // ACGT-only variant totals.
  u64 num_input_variants{0};
  u64 num_output_variants{0};
  u64 num_pass_variants{0};
  u64 num_pass_snvs{0};
  u64 num_pass_insertions{0};
  u64 num_pass_deletions{0};

  // Genotype breakdown over PASS ACGT-only variants.
  // GT=0/1
  u64 num_pass_biallelic{0};
  // GT=1/2
  u64 num_pass_multiallelic{0};
  // GT=1/1 or haploid GT=1
  u64 num_pass_homozygous{0};

  // Ti/Tv over PASS SNV ALT alleles.
  u64 num_transitions{0};
  u64 num_transversions{0};

  // Per-FILTER breakdown over non-PASS ACGT-only variants, keyed by FILTER tag.
  StrMap<u64> non_pass_filter_counts;

  /**
   * @brief Accumulate output-side counters for a single output VCF record.
   * Updates output/pass totals, the PASS variant-type and genotype breakdowns,
   * PASS Ti/Tv, and the per-FILTER breakdown. Each ACGT-only ALT allele counts as one variant.
   * @param record Output VCF record.
   */
  void AccumulateOutputRecord(const io::VcfRecord& record);

 private:
  /**
   * @brief Accumulate the PASS genotype breakdown (biallelic/multi-allelic/homozygous) for a record's ACGT variants.
   * Records without a genotype field are skipped.
   * @param record PASS output VCF record.
   * @param num_variants Number of ACGT-only ALT alleles on the record.
   */
  void AccumulatePassGenotype(const io::VcfRecord& record, u64 num_variants);

  /**
   * @brief Accumulate PASS SNV/insertion/deletion and Ti/Tv counts over the ACGT-only ALT alleles of a record.
   * Each ACGT-only ALT allele is right-trimmed and classified independently. Equal-length multi-nucleotide
   * substitutions remaining after trimming are neither SNVs nor indels.
   * @param record PASS output VCF record.
   */
  void AccumulatePassVariantTypes(const io::VcfRecord& record);

  /**
   * @brief Accumulate the per-FILTER breakdown for a non-PASS record.
   * Each ACGT-only ALT allele counts once per non-PASS FILTER tag; `PASS`/`.` tags are skipped.
   * @param record Non-PASS output VCF record.
   * @param num_variants Number of ACGT-only ALT alleles on the record.
   */
  void AccumulateNonPassFilters(const io::VcfRecord& record, u64 num_variants);

 public:
  /// Non-PASS output variants: every ACGT-only output variant that is not PASS.
  u64 NumNonPassVariants() const {
    return num_output_variants - num_pass_variants;
  }

  /// Heterozygous PASS variants: the GT=0/1 (biallelic) and GT=1/2 (multi-allelic) entries.
  u64 NumPassHeterozygous() const {
    return num_pass_biallelic + num_pass_multiallelic;
  }
};

/**
 * @brief Classify a single SNV substitution as a transition.
 * Transitions are purine<->purine (A<->G) or pyrimidine<->pyrimidine (C<->T).
 * @param ref Single-base reference allele.
 * @param alt Single-base alternate allele.
 * @return true if the substitution is a transition, false if a transversion.
 */
bool IsTransition(const std::string& ref, const std::string& alt);

/**
 * @brief Count input variants whose REF and ALT alleles contain only ACGT bases.
 * @param vcf_file Path to the input VCF.
 * @return Total number of ACGT-only ALT alleles across all records.
 */
u64 CountInputVariants(const fs::path& vcf_file);

/**
 * @brief Write VCF QC metrics to a key-value TSV file.
 * Missing/uncomputable ratios are written as `NA`. The file includes a
 * `##RocheCommandLine` metadata header and, when non-empty, a `##Sample` line.
 * @param output_file Destination TSV path (parent directories are created).
 * @param metrics Aggregated metrics.
 * @param command_line Command-line metadata for the header.
 * @param sample_name Sample name for the `##Sample` header line; for multi-sample inputs this is the
 *        comma-separated list of sample names. An empty value emits no `##Sample` line.
 * @param predicted_sex Predicted biological sex, emitted as the `predicted_sex` metric.
 * @param sample_type Final sample type after any auto-detection. `kAutoDetect` is emitted as `NA`.
 */
void WriteVcfMetricsFile(const fs::path& output_file,
                         const VcfMetrics& metrics,
                         const io::CommandLineInfo& command_line,
                         std::string_view sample_name,
                         sex_predict::Sex predicted_sex,
                         SampleType sample_type);

}  // namespace xoos::svc
