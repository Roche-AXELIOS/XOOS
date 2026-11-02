#pragma once

#include "sex-predictor-cascade.h"

namespace xoos::sex_predict {

/// Confidence for a measurement already known to clear `gate`: `kHigh` if it
/// clears by at least `kHighConfidenceMarginMultiplier`, `kMedium` otherwise.
/// See `SexConfidence`'s doc comment for the full table and why this margin
/// -- not a second calibrated threshold -- is what separates the two.
SexConfidence ConfidenceFromMargin(f64 value, f64 gate);

/// Cascade step 5: resolve from the measured heterozygosity, given chrY is
/// already known absent. Heterozygosity can CONFIRM two copies of chrX, or
/// fail to -- it deliberately cannot overturn the ratio call, because the two
/// outcomes are not symmetric: heterozygosity present is near-conclusive
/// evidence of two copies, while heterozygosity absent is consistent with an
/// LOY male, a Turner/LoX female, a tumor with chrX LOH, AND a perfectly
/// ordinary female cell line that lost an X in culture (HCC1395; see the
/// README). So an unconfirmed chrX keeps the ratio's answer and carries the
/// caveat in `reason`, where a caller that cares about LOY can act on it.
SexPrediction ResolveFromHeterozygosity(const SexMeasurements& measurements, Sex ratio_sex, SexPredictionPolicy policy);

}  // namespace xoos::sex_predict
