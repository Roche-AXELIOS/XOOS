#pragma once

#include <optional>
#include <string>

#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>

namespace xoos::cnc {
const std::string kSampleMetadataOptionsDefaultSampleId = "sample";

struct SampleMetadataOptions {
  std::string sample_id;
  sex_predict::SexSelection sex_selection = sex_predict::SexSelection::kAutoDetect;
  std::optional<sex_predict::Sex> sex;
  // Name of the normal sample in the VCF (used by VCF parsing and seg-to-vcf).
  std::optional<std::string> normal_sample_name;
  // Name of the tumor sample in the VCF (somatic T/N or tumor-only workflows).
  std::optional<std::string> tumor_sample_name;

  // Somatic Purity/Ploidy
  std::optional<f64> purity;
  std::optional<f64> ploidy;
};
}  // namespace xoos::cnc
