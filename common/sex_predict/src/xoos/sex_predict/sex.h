#pragma once

#include <concepts>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

#include <xoos/types/int.h>

namespace xoos::sex_predict {

enum class Sex : u8 {
  kUnknown,
  kMale,
  kFemale
};

/// Configures how a resolved Sex value is obtained.
enum class SexSelection : u8 {
  kAutoDetect,
  kMale,
  kFemale,
  kSkipPrediction
};

enum class SexChromHandling : u8 {
  kAllChrom,
  kNoSexChrom,
  kOnlySexChrom,
  kNoChromY
};

/// Return a capitalized human-readable name ("Male", "Female", "Unknown").
std::string GetDescriptionForSex(Sex sex);

/// Resolve a SexSelection to a Sex value.
///
/// Explicit selections map directly (kMale/kFemale) and kSkipPrediction resolves to
/// kUnknown. kAutoDetect invokes @p predict, which each module supplies to run its own
/// prediction (e.g. from coverage or a BAM). Templated so callers can pass any zero-arg
/// callable returning Sex without a std::function allocation; the predictor is taken by
/// value (the standard idiom for callables) and defined in the header as permitted for
/// template functions.
template <std::invocable Predictor>
  requires std::same_as<std::invoke_result_t<Predictor>, Sex>
Sex ResolveSex(const SexSelection selection, Predictor predict) {
  switch (selection) {
    case SexSelection::kAutoDetect:
      return predict();
    case SexSelection::kMale:
      return Sex::kMale;
    case SexSelection::kFemale:
      return Sex::kFemale;
    case SexSelection::kSkipPrediction:
      return Sex::kUnknown;
    default:
      // Unreachable: the switch above is exhaustive over SexSelection. std::runtime_error is thrown
      // directly (rather than xoos::error::Error) to keep this public header free of the fmt include.
      throw std::runtime_error("Invalid sex selection");
  }
}

/// A region string has the form "contig" or "contig:start-end" (e.g. "chrX" or
/// "chrX:1-100"). Only the contig name is compared, using genome::IsEqualContig;
/// coordinates, if present, are ignored.

/// Check if the contig of a region string is a sex chromosome (chrX, chrY, X, Y).
bool IsInAllosome(std::string_view region);

/// Check if the contig of a region string is chromosome Y (chrY, Y).
bool IsInChromY(std::string_view region);

}  // namespace xoos::sex_predict
