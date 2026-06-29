#pragma once

#include <optional>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "io/alignment.h"

namespace xoos::read_collapser {

// Base quality for two concordant reads after deconvolution.
// This value is derived from the qscore-calculator binomial model.
constexpr u8 kConcordantDeconvolvedBaseQual = 40;

/**
 * Get the discordant duplex error rate of a duplex read, which is
 * defined as the number of discordant bases divided by the total
 * duplex bases.
 *
 * @param record The pointer to the BAM record.
 *
 * @return The discordant duplex error rate (between 0 and 1) if the read is a duplex read, otherwise `std::nullopt`.
 */
std::optional<f64> GetDiscordantDuplexErrorRate(const bam1_t* record);

/**
 * Update the base qualities of an unmapped, duplex read by replacing any concordant base qualities with the deconvolved
 * quality score.
 *
 * @param read The BAM record.
 */
void UpdateUnmappedDeconvolvedQualities(const bam1_t* read);

/**
 * Get the average cluster size for a deconvolved read based on its base types.
 *
 * @param read The BAM record.
 *
 * @return The average cluster size for the deconvolved read.
 */
u32 GetUnmappedDeconvolvedAvgClusterSize(const bam1_t* read);

/**
 * Given a vector of duplex reads, deconvolve them into two separate reads based on the YC tag.
 *
 * For each read, it decodes the YC tag and generates
 * - Two separate BAM records if the read has a duplex segment
 * - One BAM record if the read is a simplex read, does not have a duplex segment, or does not have a YC tag
 *
 * @param original_alignments The vector of original alignments to deconvolve.
 * @param is_parent_parent If the reads are produced using the linearly amplified parent-parent library preparation.
 *                         Setting this correctly will ensure that the R1 and R2 reads are assigned to the correct
 * strands.
 */
vec<AlignmentPtr> DeconvolveDuplexReads(const vec<AlignmentPtr>& original_alignments, bool is_parent_parent);

}  // namespace xoos::read_collapser
