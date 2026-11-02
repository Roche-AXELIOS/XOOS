#include "gc-correct/gc-correct-algorithm.h"

#include <unordered_set>

#include <xoos/gc_correct/gc-correct.h>
#include <xoos/genome/region.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {

namespace {

// Minimum mappability for a region to be used to train the GC-correction model.
constexpr f64 kMinTrainingMappability = 1.0;

// Build the full-size gc_correct coverage record list, one entry per region, in input order.
// average_coverage is left at its default: xoos::gc_correct recomputes it internally from
// total_coverage and region size rather than reading it from the input record.
gc_correct::CoverageRecords BuildCoverageRecords(const CoverageRecords& coverage) {
  gc_correct::CoverageRecords records;
  records.reserve(coverage.region.size());
  for (size_t i = 0; i < coverage.region.size(); ++i) {
    // This caller's coverage files always give an explicit "chr:start-end" range, so restrict the
    // parser to that form to fail fast on malformed input instead of silently producing a
    // whole-chromosome span.
    const auto parsed = genome::ParseSamtoolsLikeRegion(coverage.region[i], genome::RegionForm::kRange);
    records.push_back(
        {.region = {.contig = parsed.chromosome, .start = ToUnsigned(parsed.start), .end = ToUnsigned(parsed.end)},
         .read_count = coverage.count[ToSigned(i)],
         .total_coverage = coverage.total_coverage[ToSigned(i)]});
  }
  return records;
}

// Build the per-region GC bias vector expected by xoos::gc_correct.
vec<f64> BuildGCBias(const AugmentedIntervals& augmented_intervals) {
  const Eigen::VectorXd& gc_bias = augmented_intervals.GetGCBias();
  vec<f64> out(ToUnsigned(gc_bias.size()));
  for (Eigen::Index i = 0; i < gc_bias.size(); ++i) {
    out[ToUnsigned(i)] = gc_bias[i];
  }
  return out;
}

// Indices of regions used to train the GC-correction model: autosome-only and mappability >= 1,
// matching this caller's legacy region-selection criteria for GC correction.
vec<size_t> GetTrainingIndices(const CoverageRecords& coverage, const AugmentedIntervals& augmented_intervals) {
  const Eigen::VectorXd& mappability = augmented_intervals.GetMappability();
  vec<size_t> training_indices;
  for (size_t i = 0; i < coverage.region.size(); ++i) {
    if (!sex_predict::IsInAllosome(coverage.region[i]) && mappability[ToSigned(i)] >= kMinTrainingMappability) {
      training_indices.push_back(i);
    }
  }
  return training_indices;
}

// xoos::gc_correct zero-fills any region it considers invalid for correction (non-positive
// average coverage or negative GC bias) instead of leaving it untouched. Restore the original raw
// values for those regions so they pass through unchanged, matching this caller's legacy behavior.
//
// average_coverage is recomputed from total_coverage and region size rather than read from
// `original.average_coverage`: that field is only populated when coverage is loaded from a BED file
// (CoverageRecords::LoadCoverageFile calls CalculateAverageCoverages), not when coverage comes from
// the in-pipeline BAM-derived path (BamToCorrectedCoverageTaskFlowGraph), where it stays default-sized
// and would be out of bounds here. The legacy Armadillo code had the same shape: its average-coverage
// array was always computed internally from total_coverage/region size, never read from CNC's raw
// struct.
void RestoreNonValidRegions(gc_correct::CoverageRecords& corrected,
                            const CoverageRecords& original,
                            const vec<std::ptrdiff_t>& valid_indices) {
  const std::unordered_set<std::ptrdiff_t> is_valid(valid_indices.begin(), valid_indices.end());
  for (size_t i = 0; i < corrected.size(); ++i) {
    if (is_valid.contains(ToSigned(i))) {
      continue;
    }
    corrected[i].read_count = original.count[ToSigned(i)];
    corrected[i].total_coverage = original.total_coverage[ToSigned(i)];
    const auto region_size = static_cast<f64>(corrected[i].region.end - corrected[i].region.start);
    corrected[i].average_coverage = corrected[i].total_coverage / region_size;
  }
}

// Convert a gc_correct coverage record list back into this caller's Eigen-vector layout.
GCCorrectResult ToGCCorrectResult(const gc_correct::CoverageRecords& records) {
  GCCorrectResult result;
  const auto n = ToSigned(records.size());
  result.count.resize(n);
  result.total_coverage.resize(n);
  result.average_coverage.resize(n);
  for (size_t i = 0; i < records.size(); ++i) {
    const auto idx = ToSigned(i);
    result.count[idx] = records[i].read_count;
    result.total_coverage[idx] = records[i].total_coverage;
    result.average_coverage[idx] = records[i].average_coverage;
  }
  return result;
}

}  // namespace

GCCorrectResult GCCorrect(const CoverageRecords& coverage,
                          const AugmentedIntervals& augmented_intervals,
                          const f64 first_span) {
  const gc_correct::CoverageRecords all_records = BuildCoverageRecords(coverage);
  const vec<f64> all_gc_bias = BuildGCBias(augmented_intervals);

  const vec<size_t> training_indices = GetTrainingIndices(coverage, augmented_intervals);
  gc_correct::CoverageRecords training_records;
  training_records.reserve(training_indices.size());
  vec<f64> training_gc_bias;
  training_gc_bias.reserve(training_indices.size());
  for (const size_t i : training_indices) {
    training_records.push_back(all_records[i]);
    training_gc_bias.push_back(all_gc_bias[i]);
  }

  const gc_correct::GCCorrectModel model =
      gc_correct::TrainGCCorrectModel(training_records, training_gc_bias, first_span);
  gc_correct::GCCorrectCoverage result = gc_correct::GCCorrectFromPrecomputedModel(model, all_records, all_gc_bias);
  RestoreNonValidRegions(result.corrected_coverage, coverage, result.valid_indices);
  return ToGCCorrectResult(result.corrected_coverage);
}

}  // namespace xoos::cnc
