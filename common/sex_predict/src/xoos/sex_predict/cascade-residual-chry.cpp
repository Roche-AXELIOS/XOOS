#include "cascade-residual-chry.h"

#include <xoos/log/logging.h>

#include "cascade-resolve.h"

namespace xoos::sex_predict {

std::optional<SexPrediction> TryResidualChrY(const SexMeasurements& measurements) {
  using enum SexEvidence;
  using enum SexUncertaintyReason;

  if (measurements.chr_y_depth_copies < kMinResidualChrYCopiesForMale) {
    return std::nullopt;
  }
  Logging::Debug(
      "sex_predict: chrY depth ({:.3f} copies) is below the ratio's Male threshold but above the residual gate "
      "({:.3f}); calling Male on residual chrY",
      measurements.chr_y_depth_copies,
      kMinResidualChrYCopiesForMale);
  return SexPrediction{
      .sex = Sex::kMale,
      .evidence = kResidualChrYDepth,
      .confidence = ConfidenceFromMargin(measurements.chr_y_depth_copies, kMinResidualChrYCopiesForMale),
      .reason = kNone,
      .measurements = measurements};
}

}  // namespace xoos::sex_predict
