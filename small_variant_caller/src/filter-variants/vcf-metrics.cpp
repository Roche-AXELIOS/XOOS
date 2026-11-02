#include "vcf-metrics.h"

#include <optional>
#include <string>
#include <string_view>

#include <fmt/format.h>

#include <xoos/enum/enum-util.h>
#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/constants.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/vec.h>
#include <xoos/util/missing-value.h>

#include "core/genotype.h"
#include "util/locked-tsv-writer.h"
#include "util/seq-util.h"

namespace xoos::svc {

namespace {

/**
 * @brief Number of ACGT-only ALT alleles in a record (each counts as one variant).
 * @param record VCF record to inspect.
 * @return Count of ALT alleles for which both REF and ALT contain only ACGT bases.
 */
u64 NumAcgtVariants(const io::VcfRecord& record) {
  if (!ContainsOnlyACTG(record.Allele(0))) {
    return 0;
  }

  const s32 num_alleles = record.NumAlleles();
  u64 count = 0;
  for (s32 i = 1; i < num_alleles; ++i) {
    if (ContainsOnlyACTG(record.Allele(i))) {
      ++count;
    }
  }
  return count;
}

/**
 * @brief Compute a ratio as an optional; empty when the denominator is zero.
 */
std::optional<f64> SafeRatio(const u64 numerator, const u64 denominator) {
  if (denominator == 0) {
    return std::nullopt;
  }
  return static_cast<f64>(numerator) / static_cast<f64>(denominator);
}

/**
 * @brief Compute a percentage as an optional; empty when the base is zero.
 */
std::optional<f64> SafePercentage(const u64 count, const u64 base) {
  if (base == 0) {
    return std::nullopt;
  }
  return 100.0 * static_cast<f64>(count) / static_cast<f64>(base);
}

/// The three value columns of a key-value metrics row. `percentage`/`denominator` default to `NA`.
struct MetricRow {
  std::string value;
  std::string percentage{util::kNA};
  std::string denominator{util::kNA};
};

/// Build one key-value row as the four output columns: metric_name, value, percentage, denominator.
vec<std::string> MakeRow(const std::string_view metric_name, const MetricRow& row) {
  return {std::string{metric_name}, row.value, row.percentage, row.denominator};
}

/// Build a count row whose percentage is relative to a named base count.
vec<std::string> MakeCountRow(const std::string_view metric_name,
                              const u64 count,
                              const u64 base,
                              const std::string_view denominator) {
  return MakeRow(metric_name,
                 {.value = fmt::format("{}", count),
                  .percentage = util::FormatOrNA(SafePercentage(count, base), "{:.2f}"),
                  .denominator = std::string{denominator}});
}

}  // namespace

bool IsTransition(const std::string& ref, const std::string& alt) {
  if (ref.size() != 1 || alt.size() != 1) {
    return false;
  }
  // A transition is an interchange within a base class: the purines A<->G or the pyrimidines C<->T.
  // Membership is case-insensitive; any non-ACGT base belongs to neither class and is not a transition.
  constexpr std::string_view kPurines = "AGag";
  constexpr std::string_view kPyrimidines = "CTct";
  const char ref_base = ref[0];
  const char alt_base = alt[0];
  if (ref_base == alt_base) {
    return false;
  }
  const bool both_purines = kPurines.contains(ref_base) && kPurines.contains(alt_base);
  const bool both_pyrimidines = kPyrimidines.contains(ref_base) && kPyrimidines.contains(alt_base);
  return both_purines || both_pyrimidines;
}

void VcfMetrics::AccumulatePassGenotype(const io::VcfRecord& record, const u64 num_variants) {
  // Genotype breakdown is variant-level, matching num_pass_variants. Each ACGT-only ALT allele on the record is
  // assigned to the record's genotype category.
  // Guard the GT read so records without a genotype field are simply skipped.
  if (!record.HasFormatFieldNoCheck(io::kGT)) {
    return;
  }
  using enum Genotype;
  switch (StringToGenotype(record.GetGTField())) {
    case kGT01:
      num_pass_biallelic += num_variants;
      break;
    case kGT12:
      num_pass_multiallelic += num_variants;
      break;
    case kGT11:
    case kGT1:
      num_pass_homozygous += num_variants;
      break;
    default:
      break;
  }
}

void VcfMetrics::AccumulatePassVariantTypes(const io::VcfRecord& record) {
  const std::string ref = record.Allele(0);
  if (!ContainsOnlyACTG(ref)) {
    return;
  }

  const s32 num_alleles = record.NumAlleles();
  for (s32 i = 1; i < num_alleles; ++i) {
    const std::string alt = record.Allele(i);
    if (!ContainsOnlyACTG(alt)) {
      continue;
    }
    const auto& [trimmed_ref, trimmed_alt] = TrimVariant(ref, alt);
    if (trimmed_ref.size() < trimmed_alt.size()) {
      ++num_pass_insertions;
      continue;
    }
    if (trimmed_ref.size() > trimmed_alt.size()) {
      ++num_pass_deletions;
      continue;
    }
    if (trimmed_ref.size() != 1) {
      continue;
    }

    ++num_pass_snvs;
    if (IsTransition(trimmed_ref, trimmed_alt)) {
      ++num_transitions;
    } else {
      ++num_transversions;
    }
  }
}

void VcfMetrics::AccumulateNonPassFilters(const io::VcfRecord& record, const u64 num_variants) {
  // Count each FILTER tag that is not `PASS`/`.`, so the breakdown holds one entry per genuine
  // failure reason even on the (malformed but valid) records that mix `PASS` with other tags.
  for (const auto& filter : record.GetFilters()) {
    if (filter == io::kFilterPASS || filter == io::kFilterDOT) {
      continue;
    }
    non_pass_filter_counts[filter] += num_variants;
  }
}

void VcfMetrics::AccumulateOutputRecord(const io::VcfRecord& record) {
  const u64 num_variants = NumAcgtVariants(record);
  num_output_variants += num_variants;

  if (num_variants == 0) {
    return;
  }
  if (!record.IsPass()) {
    AccumulateNonPassFilters(record, num_variants);
    return;
  }

  num_pass_variants += num_variants;
  AccumulatePassGenotype(record, num_variants);
  AccumulatePassVariantTypes(record);
}

u64 CountInputVariants(const fs::path& vcf_file) {
  const io::VcfReader reader(vcf_file);
  u64 count = 0;
  while (const auto record = reader.GetNextRecord()) {
    count += NumAcgtVariants(*record);
  }
  return count;
}

void WriteVcfMetricsFile(const fs::path& output_file,
                         const VcfMetrics& metrics,
                         const io::CommandLineInfo& command_line,
                         const std::string_view sample_name,
                         const sex_predict::Sex predicted_sex,
                         const SampleType sample_type) {
  // LockedTsvWriter creates the parent directory on open and throws error::Error on failure.
  LockedTsvWriter writer(output_file);

  // Metadata header: the ##RocheCommandLine and ##Sample lines are guarded independently to preserve
  // output for internal callers that supply no command line (empty version) and/or no sample name.
  // Metadata must be appended before any rows: the writer flushes buffered rows and writes metadata
  // comment lines directly to the underlying stream (bypassing CSV quoting).
  if (!command_line.version.empty()) {
    writer.AppendMetadata(command_line);
  }
  if (!sample_name.empty()) {
    io::ReadGroup read_group;
    read_group.sample_name = sample_name;
    writer.AppendSampleMetadata(read_group);
  }

  // QC pass/fail: fail when there are no output variants or no PASS variants.
  const bool qc_pass = metrics.num_output_variants > 0 && metrics.num_pass_variants > 0;

  // Writer-bound emitters so each metric line reads as a single builder call.
  const auto add_row = [&writer](const std::string_view metric_name, const MetricRow& value) {
    writer.AppendRow(MakeRow(metric_name, value));
  };
  const auto add_count_row =
      [&writer](
          const std::string_view metric_name, const u64 count, const u64 base, const std::string_view denominator) {
        writer.AppendRow(MakeCountRow(metric_name, count, base, denominator));
      };

  // Column header row (snake_case, unprefixed).
  writer.AppendRow({std::string{kMetricColumnName},
                    std::string{kMetricColumnValue},
                    std::string{kMetricColumnPercentage},
                    std::string{kMetricColumnDenominator}});

  // Overall totals.
  add_row(kMetricQcStatus, {.value = std::string{qc_pass ? kQcStatusPass : kQcStatusFail}});
  add_row(kMetricPredictedSex, {.value = enum_util::FormatEnumName(predicted_sex)});
  add_row(kMetricSampleType,
          {.value = sample_type == SampleType::kAutoDetect ? std::string{util::kNA}
                                                           : enum_util::FormatEnumName(sample_type)});
  add_row(kMetricNumInputVariants, {.value = fmt::format("{}", metrics.num_input_variants)});
  add_count_row(
      kMetricNumOutputVariants, metrics.num_output_variants, metrics.num_input_variants, kMetricNumInputVariants);
  add_count_row(
      kMetricNumPassVariants, metrics.num_pass_variants, metrics.num_output_variants, kMetricNumOutputVariants);
  add_count_row(
      kMetricNumNonPassVariants, metrics.NumNonPassVariants(), metrics.num_output_variants, kMetricNumOutputVariants);
  add_count_row(kMetricNumPassSnvs, metrics.num_pass_snvs, metrics.num_pass_variants, kMetricNumPassVariants);
  add_count_row(
      kMetricNumPassInsertions, metrics.num_pass_insertions, metrics.num_pass_variants, kMetricNumPassVariants);
  add_count_row(kMetricNumPassDeletions, metrics.num_pass_deletions, metrics.num_pass_variants, kMetricNumPassVariants);

  // PASS genotype breakdown (variant-level), relative to PASS variants.
  add_count_row(kMetricNumPassBiallelic, metrics.num_pass_biallelic, metrics.num_pass_variants, kMetricNumPassVariants);
  add_count_row(
      kMetricNumPassMultiallelic, metrics.num_pass_multiallelic, metrics.num_pass_variants, kMetricNumPassVariants);
  add_count_row(
      kMetricNumPassHomozygous, metrics.num_pass_homozygous, metrics.num_pass_variants, kMetricNumPassVariants);
  add_count_row(
      kMetricNumPassHeterozygous, metrics.NumPassHeterozygous(), metrics.num_pass_variants, kMetricNumPassVariants);

  // Ratios.
  add_row(kMetricHetHomRatio,
          {.value = util::FormatOrNA(SafeRatio(metrics.NumPassHeterozygous(), metrics.num_pass_homozygous), "{:.2f}")});
  add_count_row(kMetricNumTransitions, metrics.num_transitions, metrics.num_pass_snvs, kMetricNumPassSnvs);
  add_count_row(kMetricNumTransversions, metrics.num_transversions, metrics.num_pass_snvs, kMetricNumPassSnvs);
  add_row(kMetricTiTvRatio,
          {.value = util::FormatOrNA(SafeRatio(metrics.num_transitions, metrics.num_transversions), "{:.2f}")});

  // Per-FILTER breakdown over non-PASS variants (std::map keeps tags sorted).
  for (const auto& [filter, count] : metrics.non_pass_filter_counts) {
    add_count_row(fmt::format("{}{}", kMetricFilterPrefix, filter),
                  count,
                  metrics.NumNonPassVariants(),
                  kMetricNumNonPassVariants);
  }
}

}  // namespace xoos::svc
