#pragma once

#include <xoos/types/int.h>

namespace xoos::demux::scoring {

/**
 * @namespace scoring
 * @brief Sequencer-specific log-odds alignment scoring constants. As opposed to edit distance, methods that use these
 * constants will using a probabilistic alignment model which should be more sensitive/specific, though if misused
 * possibly slower.
 *
 * These values model the error profile of our sequencing technology as integer log-odds ratios.
 * They are additive (log-space) and can be used by any alignment or classification engine (DFA classifier,
 * Smith-Waterman, etc.) that needs to score matches, mismatches, insertions, and deletions against a reference.
 *
 * The values are currently best guess hand-tuned based on empirical observations of SBX data.
 * They are trainable: We can use expected adapter elements to calibrate directly
 * but we keep them as integers to constrain the delta range for speed/size in engines like the DFA.
 */

// Match score of 2 (in bits) mathematically reflects that a true biological match is 4 times more likely to produce
// identical bases than random chance (1.0 (match) vs 0.25 (random) probability), because log2(1.0 / 0.25) = 2.
// The log odds here is:
// Likelihood that matching base exists if a real sequence match = 100%,
// divided by,
// Likelihood that the matching base exists if part of random sequence = 1/4 = 25%
// 1/(1/4) = 4 taking the log2(4) of this is 2
inline constexpr s32 kMatch = 2;

// Penalties are scores that represent the likelihood in the same way as the match score, but they are negative
// because they are log odds ratios less than 1. They basically decrease probability we think it is an overall match.
// For example, the log odds ratio for substitutions:
// Likelihood that base isn't matching because of a mismatch due to sequencer noise = x,
// divided by,
// Likelihood that non-matching base exists if part of a random sequence = 3/4 = 0.75
// log2(x/0.75) = -4
// Solving here for x = (2^-4)*(0.75) = 0.046875 (or 1 in ~21 chance of mismatch due to sequencer noise)
inline constexpr s32 kSubstitution = -4;
// (2^-3)*0.75 = ~0.09375 (1 in ~11)
inline constexpr s32 kInsertion = -3;
// (2^-5)*0.75 = ~0.0234375 (1 in ~43)
inline constexpr s32 kDeletion = -5;

// TODO: Reference-guided homopolymer insertion scoring. If a reference position is within a homopolymer run
//       (>= threshold consecutive identical bases) we cane soften the insertion penalty for that base.
//       This reflects that our basecaller may be more likely to produce insertion errors at reference homopolymer
//       sites. This will not address de novo (not caused by existing homopolymer) homopolymer insertion expansions
//       that may also be occurring. In our DFA, the change will be added to the AdvanceNfa function. Also de novo
//       homopolymer events can be added to the DFA, but will probably increase size more than this proposed change.
// inline constexpr s32 kHomopolymerLengthThreshold = 3;
// inline constexpr s32 kHomopolymerInsertion = -1;
// Some useful aliases
// The minimum possible single-step penalty (most negative value in the scoring scheme).
inline constexpr s32 kMinPenalty = kDeletion;

}  // namespace xoos::demux::scoring
