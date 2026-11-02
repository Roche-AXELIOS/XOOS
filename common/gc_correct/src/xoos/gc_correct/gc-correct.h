#pragma once

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "xoos/gc_correct/coverage-record.h"

namespace xoos::gc_correct {

const f64 kGCCorrectSecondLoessSpan = 0.3;
const s32 kGCCorrectLoessDegree = 2;
const f64 kGCCorrectLoessGCMin = 0;
const f64 kGCCorrectLoessGCMax = 1;
const f64 kGCCorrectLoessGCStep = 0.001;

// Result struct containing GC-corrected coverage data, scale factors, and indices of valid data points used in
// correction.
struct GCCorrectCoverage {
  vec<CoverageRecord> corrected_coverage;
  vec<f64> gc_bins;
  vec<f64> smoothed_coverage_at_gc_bins;
  vec<f64> scale_factor;
  vec<std::ptrdiff_t> valid_indices;
};

struct GCCorrectInput {
  vec<f64> counts;
  vec<f64> total_coverage;
  vec<f64> average_coverage;
  vec<f64> region_sizes;
};

// Intermediate data structure containing filtered regions and their metrics for LOESS model training and application.
struct GCCorrectLoessInput {
  // Indices of regions with valid coverage (>0) and GC bias (>=0)
  vec<std::ptrdiff_t> is_valid_idxs;
  // Indices of valid regions without coverage and GC bias outliers
  vec<std::ptrdiff_t> is_ideal_idxs;
  // Average coverage values for ideal regions used to train LOESS model
  vec<f64> ideal_average_coverage;
  // GC bias values for ideal regions used to train LOESS model
  vec<f64> ideal_gc_bias;
  // Average coverage values for all valid regions (used in LOESS smoothing)
  vec<f64> valid_average_coverage;
  // GC bias values for all valid regions (used in LOESS smoothing)
  vec<f64> valid_gc_bias;
};

struct GCCorrectModel {
  vec<f64> x;
  vec<f64> y;
  f64 mean;
  GCCorrectInput input;
  GCCorrectLoessInput loess_input;
};

/**
 * @brief Train a GC correction model and apply it on the same data to produce GC-corrected coverage
 * @param coverage - vector of CoverageRecord objects containing coverage data for genomic regions.
 * @param gc_bias  - vector of GC bias values corresponding to the genomic regions described by coverage.
 * @param first_span - span parameter > 0, <= 1 for the first round of LOESS
 * smoothing. Corresponds to how the fraction data used to smooth each point.
 * Lower values are faster but might be more inaccurate. Higher values might
 * result in over-fitting. For small datasets use 1.
 * @return GCCorrectCoverage struct containing GC-corrected coverage data, scale factors, and indices of valid data
 * points used in correction.
 */
GCCorrectCoverage GCCorrect(const vec<CoverageRecord>& coverage, const vec<f64>& gc_bias, f64 first_span);

/**
 * @brief Train a GC correction model using LOESS smoothing. No correction is applied.
 * @param coverage - vector of CoverageRecord objects containing coverage data for genomic regions.
 * @param gc_bias  - vector of GC bias values corresponding to the genomic regions described by coverage.
 * @param first_span - span parameter > 0, <= 1 for the first round of LOESS
 * smoothing. Corresponds to how the fraction data used to smooth each point.
 * Lower values are faster but might be more inaccurate. Higher values might
 * result in over-fitting. For small datasets use 1.
 * @return GCCorrectModel struct containing the trained LOESS model parameters and relevant intermediate data required
 * for generating the model.
 */
GCCorrectModel TrainGCCorrectModel(const vec<CoverageRecord>& coverage, const vec<f64>& gc_bias, f64 first_span);

/**
 * @brief Apply a pre-computed GC corrected model to new data
 * @param model     - GCCorrectModel struct containing the pre-computed GC correction model. Generate this using
 * TrainGCCorrectModel
 * @param coverage  - vector of CoverageRecord objects containing coverage data for genomic regions.
 * @param gc_bias   - vector of GC bias values corresponding to the genomic regions described by coverage.
 * @return GCCorrectCoverage struct containing GC-corrected coverage data
 */
GCCorrectCoverage GCCorrectFromPrecomputedModel(const GCCorrectModel& model,
                                                const vec<CoverageRecord>& coverage,
                                                const vec<f64>& gc_bias);

}  // namespace xoos::gc_correct
