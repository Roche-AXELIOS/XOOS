#pragma once

#include <optional>

#include "sex-predictor-cascade.h"

namespace xoos::sex_predict {

/// Cascade step 3b: residual chrY depth, evidence of partial (not complete)
/// loss-of-Y. Returns a Male call when `measurements.chr_y_depth_copies`
/// clears `kMinResidualChrYCopiesForMale`, or `std::nullopt` when it doesn't
/// (the cascade continues to the next step). See the README's "Residual
/// chrY" section for the full derivation.
std::optional<SexPrediction> TryResidualChrY(const SexMeasurements& measurements);

}  // namespace xoos::sex_predict
