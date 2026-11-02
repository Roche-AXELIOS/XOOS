#pragma once

#include <string>

#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::svc::phaser {

/** @brief A variant record used by the phaser. */
struct Variant {
  std::string chrom;
  // 0-based
  s64 pos = 0;
  std::string ref;
  std::string alt;
  std::string filter;
  // index in the original VCF for output mapping
  s32 vcf_record_index = -1;
  // true if FILTER is PASS or "." (unfiltered)
  bool is_pass = true;

  /** @brief Return whether the variant is a single-nucleotide variant. */
  bool IsSnv() const {
    return ref.size() == 1 && alt.size() == 1;
  }

  /** @brief Return whether the variant is a multi-nucleotide variant. */
  bool IsMnv() const {
    return ref.size() == alt.size() && ref.size() > 1;
  }

  /** @brief Return whether the alternate allele is longer than the reference allele. */
  bool IsInsertion() const {
    return alt.size() > ref.size();
  }

  /** @brief Return whether the reference allele is longer than the alternate allele. */
  bool IsDeletion() const {
    return ref.size() > alt.size();
  }

  /**
   * @brief Return the exclusive reference-coordinate end of the variant.
   * @return Zero-based end position, computed as `pos + ref.size()`.
   */
  s64 EndPos() const {  // NOSONAR(cpp:S5536) — used throughout phase-engine.cpp and vcf-io.cpp
    return pos + static_cast<s64>(ref.size());
  }
};

/** @brief A group of nearby variants that are candidates for phasing. */
struct PhaseGroup {
  // Indices into the variant list.
  vec<s32> variant_indices;
  std::string chrom;
};

/** @brief Result of read-based phasing between two variants. */
struct PhasingResult {
  s32 variant_idx_a = -1;
  s32 variant_idx_b = -1;
  // reads with both ALTs
  s32 cis_count = 0;
  // reads with one ALT and one REF
  s32 trans_count = 0;
  // reads with both REFs
  s32 ref_count = 0;
  s32 total_spanning = 0;
  bool is_phased = false;
  // true = cis (both ALTs on same haplotype), false = trans.
  bool is_cis = true;

  /**
   * @brief Return the fraction of informative reads supporting cis.
   * @return Cis reads divided by cis plus trans reads, or zero when none are informative.
   */
  f32 CisFraction() const {
    const s32 informative = cis_count + trans_count;
    return informative > 0 ? static_cast<f32>(cis_count) / static_cast<f32>(informative) : 0.0F;
  }

  /**
   * @brief Return the fraction of informative reads supporting trans.
   * @return Trans reads divided by cis plus trans reads, or zero when none are informative.
   */
  f32 TransFraction() const {
    const s32 informative = cis_count + trans_count;
    return informative > 0 ? static_cast<f32>(trans_count) / static_cast<f32>(informative) : 0.0F;
  }
};

/** @brief A connected set of phased variants sharing a phase set. */
struct PhaseBlock {
  vec<s32> variant_indices;
  // Haplotype assignment for each variant (parallel to variant_indices).
  // 0 = haplotype A (PGT 0|1), 1 = haplotype B (PGT 1|0).
  // Empty when haplotype assignment is not available.
  vec<s32> haplotypes;
};

}  // namespace xoos::svc::phaser
