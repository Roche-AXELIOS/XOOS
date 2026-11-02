#pragma once

#include <string>
#include <utility>

#include <htslib/sam.h>

#include <xoos/io/alignment-reader.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "allele.h"
#include "params.h"
#include "variant.h"

namespace xoos::svc::phaser {

/**
 * @brief Group nearby variants into candidate phase groups based on distance.
 * @param variants Variants to group, in any order.
 * @param max_distance Maximum allowed distance between consecutive variants in a group.
 * @return Candidate groups containing at least two variants.
 */
vec<PhaseGroup> BuildPhaseGroups(const vec<Variant>& variants, s32 max_distance);

/**
 * @brief Determine the allele observed at a variant position on a BAM read.
 * @param seq BAM-encoded read sequence.
 * @param cigar BAM CIGAR operations for the read.
 * @param n_cigar Number of CIGAR operations in `cigar`.
 * @param read_pos Reference position where the read alignment starts.
 * @param var Variant whose allele should be determined.
 * @return The observed reference, alternate, or unknown allele.
 */
Allele DetermineAllele(const u8* seq, const u32* cigar, s32 n_cigar, s64 read_pos, const Variant& var);

/**
 * @brief Check whether two variants overlap in reference coordinates.
 * @param a First variant.
 * @param b Second variant.
 * @return `true` when both variants are on the same chromosome and their reference spans overlap.
 */
bool VariantsOverlap(const Variant& a, const Variant& b);

/**
 * @brief Determine alleles for an overlapping variant pair using haplotype sequences.
 * @param seq BAM-encoded read sequence.
 * @param cigar BAM CIGAR operations for the read.
 * @param n_cigar Number of CIGAR operations in `cigar`.
 * @param read_pos Reference position where the read alignment starts.
 * @param var_a First variant in the pair.
 * @param var_b Second variant in the pair.
 * @return The alleles observed for `var_a` and `var_b`, respectively.
 */
std::pair<Allele, Allele> DetermineAllelesOverlapping(
    const u8* seq, const u32* cigar, s32 n_cigar, s64 read_pos, const Variant& var_a, const Variant& var_b);

/** @brief Thresholds controlling read-based phasing decisions for a variant pair. */
struct PhasingThresholds {
  s32 min_shared_reads = 2;
  s32 min_cis_reads = 2;
  f32 min_cis_fraction = kDefaultMinCisFraction;
  s32 min_base_quality = 0;
};

/**
 * @brief Perform read-based phasing for a pair of variants using a pre-opened alignment reader.
 * @param bam_handle Open BAM reader containing the reads to inspect.
 * @param tumor_rg Tumor read-group identifier whose reads should be counted.
 * @param var_a First variant in the pair.
 * @param var_b Second variant in the pair.
 * @param thresholds Read-count, fraction, and base-quality thresholds.
 * @return Read evidence and the resulting cis/trans phasing decision.
 */
PhasingResult PhaseVariantPair(const io::AlignmentReader& bam_handle,
                               const std::string& tumor_rg,
                               const Variant& var_a,
                               const Variant& var_b,
                               const PhasingThresholds& thresholds);

/**
 * @brief Run single-threaded phasing on all candidate groups.
 * @param bam_path Indexed BAM containing tumor reads.
 * @param tumor_rg Tumor read-group identifier whose reads should be counted.
 * @param variants Variants referenced by `groups`.
 * @param groups Candidate variant groups to phase.
 * @param params Phasing thresholds and output-related parameters.
 * @return Connected phase blocks derived from the pairwise results.
 */
vec<PhaseBlock> RunPhasing(const fs::path& bam_path,
                           const std::string& tumor_rg,
                           const vec<Variant>& variants,
                           const vec<PhaseGroup>& groups,
                           const Params& params);

/**
 * @brief Run phasing on candidate groups using multiple worker threads.
 * @param bam_path Indexed BAM containing tumor reads.
 * @param tumor_rg Tumor read-group identifier whose reads should be counted.
 * @param variants Variants referenced by `groups`.
 * @param groups Candidate variant groups to phase.
 * @param params Phasing thresholds, thread count, and output-related parameters.
 * @return Connected phase blocks derived from the pairwise results.
 */
vec<PhaseBlock> RunPhasingParallel(const fs::path& bam_path,
                                   const std::string& tumor_rg,
                                   const vec<Variant>& variants,
                                   const vec<PhaseGroup>& groups,
                                   const Params& params);

}  // namespace xoos::svc::phaser
