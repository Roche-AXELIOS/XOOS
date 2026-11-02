#pragma once

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "core/variant-info.h"

namespace xoos::svc {

/**
 * @brief Add the germline-mitochondria AF/depth FILTER header lines (`low_af`, `low_dp`).
 * @param hdr VCF header to modify
 * @param min_af Minimum allele frequency threshold (for the description)
 * @param min_dp Minimum depth threshold (for the description)
 */
void AddMitochondriaFilterHeaderLines(const io::VcfHeaderPtr& hdr, f32 min_af, u32 min_dp);

/**
 * @brief Add the optional germline-mitochondria `blocklist` FILTER header line.
 * @param hdr VCF header to modify
 */
void AddMitochondriaBlocklistFilterHeaderLine(const io::VcfHeaderPtr& hdr);

/**
 * @brief Add the FORMAT header lines for the data fields written by AF/depth filtering.
 *
 * Adds AD (allele depth) and DP (read depth) lines, plus the duplex count lines (ADC/ADS/ADD/ADL)
 * for duplex protocols. These mirror the data FORMAT fields produced by the germline workflow.
 * htslib ignores duplicate IDs, so it is safe to call on a header cloned from the input VCF.
 *
 * @param hdr VCF header to modify
 * @param is_duplex_protocol Whether the sequencing protocol is duplex (controls duplex count lines)
 */
void AddMitochondriaDataFormatLines(const io::VcfHeaderPtr& hdr, bool is_duplex_protocol);

/**
 * @brief Apply allele-frequency and depth filters to a mitochondria variant record.
 *
 * Sets the genotype based on `duplex_af` (`1/1` if `duplex_af >= min_hom_alt_af`, otherwise `0/1`),
 * populates the AD/DP (and, for duplex protocols, the duplex count) FORMAT fields, and assigns the
 * `low_af` and/or `low_dp` FILTER tags when the variant's `duplex_af` / `duplex_dp` fall below the
 * provided thresholds. A passing variant is set to `PASS`.
 *
 * @param record VCF record to be updated with filter results, genotype, and data FORMAT fields
 * @param var_feat Variant allele BAM features (provides `duplex_af`)
 * @param ref_feat Reference allele BAM features (provides `duplex_dp`)
 * @param min_af Minimum allele frequency threshold
 * @param min_dp Minimum depth threshold
 * @param min_hom_alt_af AF threshold at or above which the genotype is homozygous ALT (`1/1`)
 * @param is_duplex_protocol Whether the sequencing protocol is duplex (controls duplex count fields)
 */
void ApplyMitochondriaFiltersToRecord(const io::VcfRecordPtr& record,
                                      const UnifiedVariantFeature& var_feat,
                                      const UnifiedReferenceFeature& ref_feat,
                                      f32 min_af,
                                      u32 min_dp,
                                      f64 min_hom_alt_af,
                                      bool is_duplex_protocol);

/**
 * @brief Populate a mitochondria record that has no BAM features with default values and tag it.
 *
 * Writes the same FORMAT field set as a normal record (GT/AD/DP/AF, plus the duplex count fields for
 * duplex protocols) but with default values — a no-call genotype (`./.`) and zero-valued depth and
 * allele-frequency fields — and sets the `missing_feature` FILTER tag. This keeps the output VCF's
 * FORMAT columns consistent across records instead of dropping variants that lack features.
 *
 * @param record VCF record to populate with default FORMAT values and the `missing_feature` FILTER
 * @param is_duplex_protocol Whether the sequencing protocol is duplex (controls duplex count fields)
 */
void SetMissingFeatureRecord(const io::VcfRecordPtr& record, bool is_duplex_protocol);

/**
 * @brief Add the mitochondria `blocklist` FILTER while preserving existing failure reasons.
 *
 * Replaces PASS when the record currently passes and otherwise appends `blocklist` to the
 * existing FILTER values.
 *
 * @param record VCF record to update
 */
void AddMitochondriaBlocklistFilter(const io::VcfRecordPtr& record);

}  // namespace xoos::svc
