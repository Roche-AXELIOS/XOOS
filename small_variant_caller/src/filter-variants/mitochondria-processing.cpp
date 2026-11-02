#include "mitochondria-processing.h"

#include <string>

#include <fmt/format.h>

#include <xoos/types/vec.h>

#include "core/filtering.h"
#include "core/genotype.h"
#include "core/vcf-fields.h"
#include "util/vcf-util.h"

namespace xoos::svc {

namespace {
constexpr auto* kFilterPass = "PASS";
}  // namespace

void AddMitochondriaFilterHeaderLines(const io::VcfHeaderPtr& hdr, const f32 min_af, const u32 min_dp) {
  hdr->AddFilterLine({.id = kFilterPass, .description = "All filters passed"});

  const auto low_af_desc = fmt::format("Allele frequency is below threshold (duplex_af < {})", min_af);
  hdr->AddFilterLine({.id = kFilterMitoLowAf, .description = low_af_desc});

  const auto low_dp_desc = fmt::format("Depth is below threshold (duplex_dp < {})", min_dp);
  hdr->AddFilterLine({.id = kFilterMitoLowDp, .description = low_dp_desc});

  // Applied when a variant has no BAM features to evaluate; the record is emitted with default
  // (zero / no-call) FORMAT values rather than dropped.
  hdr->AddFilterLine({.id = kFilteringMissingFeatureId, .description = kFilteringMissingFeatureDesc});
}

void AddMitochondriaBlocklistFilterHeaderLine(const io::VcfHeaderPtr& hdr) {
  hdr->AddFilterLine({.id = kFilterMitoBlocklist, .description = kFilterMitoBlocklistDesc});
}

void AddMitochondriaDataFormatLines(const io::VcfHeaderPtr& hdr, const bool is_duplex_protocol) {
  using enum io::FieldType;
  hdr->AddFormatLine({.id = kFieldGt, .description = "Genotype", .number = io::kNumberOne, .type = kString});
  // AD carries one value per allele (REF + ALT), i.e. VCF Number=R.
  hdr->AddFormatLine({.id = kAlleleDepthId,
                      .description = "Allelic depths for the ref and alt alleles in the order listed",
                      .number = io::kNumberR,
                      .type = kInteger});
  hdr->AddFormatLine({.id = kFieldDp,
                      .description = "Approximate read depth after filtering",
                      .number = io::kNumberOne,
                      .type = kInteger});
  // AF carries one value per ALT allele, i.e. VCF Number=A.
  hdr->AddFormatLine(
      {.id = kFieldAf, .description = "Allele frequency", .number = io::kNumberEachAllele, .type = kFloat});
  if (is_duplex_protocol) {
    hdr->AddFormatLine({.id = kDuplexConcordantCountsId,
                        .description = kDuplexConcordantCountsDesc,
                        .number = io::kNumberR,
                        .type = kInteger});
    hdr->AddFormatLine({.id = kDuplexSimplexCountsId,
                        .description = kDuplexSimplexCountsDesc,
                        .number = io::kNumberR,
                        .type = kInteger});
    hdr->AddFormatLine({.id = kDuplexDiscordantCountsId,
                        .description = kDuplexDiscordantCountsDesc,
                        .number = io::kNumberR,
                        .type = kInteger});
    hdr->AddFormatLine(
        {.id = kDuplexLowbqCountsId, .description = kDuplexLowbqCountsDesc, .number = io::kNumberR, .type = kInteger});
  }
}

void ApplyMitochondriaFiltersToRecord(const io::VcfRecordPtr& record,
                                      const UnifiedVariantFeature& var_feat,
                                      const UnifiedReferenceFeature& ref_feat,
                                      const f32 min_af,
                                      const u32 min_dp,
                                      const f64 min_hom_alt_af,
                                      const bool is_duplex_protocol) {
  // `duplex_dp` is shared by all variants at a position and lives on the reference feature.
  const f64 duplex_dp = ref_feat.duplex_dp;
  const f64 duplex_af = var_feat.duplex_af;

  const bool pass_af = duplex_af >= min_af;
  const bool pass_dp = duplex_dp >= min_dp;

  // Set genotype based on duplex_af: 1/1 if >= min_hom_alt_af, otherwise 0/1
  static const auto kGT11String = GenotypeToString(Genotype::kGT11);
  static const auto kGT01String = GenotypeToString(Genotype::kGT01);
  const std::string gt = (duplex_af >= min_hom_alt_af) ? kGT11String : kGT01String;
  record->SetGTField(gt);

  // Populate the AD/DP (and duplex count) FORMAT fields shared with the germline workflow.
  SetRecordDataFormatFields(record, ref_feat, var_feat, is_duplex_protocol);

  // Write the computed allele frequency as a per-sample FORMAT field (one value for the single
  // ALT allele). This is mitochondria-specific: the shared germline data-field writer does not
  // emit FORMAT/AF, so the germline output is left unchanged.
  record->SetFormatField(kFieldAf, vec<f32>{static_cast<f32>(duplex_af)});

  if (pass_af && pass_dp) {
    record->SetFilter(kFilterPass);
    return;
  }

  // Clear existing filters and add failed filters
  record->ClearFilters();
  if (!pass_af) {
    record->AddFilter(kFilterMitoLowAf);
  }
  if (!pass_dp) {
    record->AddFilter(kFilterMitoLowDp);
  }
}

void SetMissingFeatureRecord(const io::VcfRecordPtr& record, const bool is_duplex_protocol) {
  // No BAM features are available for this variant, so write the same FORMAT field set as a normal
  // record but with default values: no-call genotype and zero-valued depth/allele-frequency fields.
  static const auto kGTNAString = GenotypeToString(Genotype::kGTNA);
  record->SetGTField(kGTNAString);

  // AD/DP (and, for duplex protocols, the duplex count fields) are zero-filled via the zero feature
  // structs shared with the germline workflow.
  SetRecordDataFormatFields(record, kZeroUnifiedReferenceFeature, kZeroUnifiedVariantFeature, is_duplex_protocol);

  // AF mirrors the mitochondria-specific write in ApplyMitochondriaFiltersToRecord (one value per ALT).
  record->SetFormatField(kFieldAf, vec<f32>{0.0F});

  record->SetFilter(kFilteringMissingFeatureId);
}

void AddMitochondriaBlocklistFilter(const io::VcfRecordPtr& record) {
  if (record->IsPass()) {
    record->SetFilter(kFilterMitoBlocklist);
    return;
  }
  record->AddFilter(kFilterMitoBlocklist);
}

}  // namespace xoos::svc
