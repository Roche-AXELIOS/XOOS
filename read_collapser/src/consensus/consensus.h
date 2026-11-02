#pragma once

#include <xoos/types/float.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "clustering/clustering.h"
#include "consensus/consensus-matrix.h"
#include "consensus/consensus-result.h"
#include "core/read-collapser-options.h"
#include "core/region-lookup.h"
#include "core/region.h"
#include "core/thread-context.h"
#include "metrics/metrics.h"

namespace xoos::read_collapser {

// MSA options.
constexpr s32 kMatchScore = 10;
constexpr s32 kMismatchPenalty = 8;
constexpr s32 kGapOpenPenalty = 80;
constexpr s32 kGapExtPenalty = 20;

// Chromosome symbol for unmapped reads.
constexpr std::string kUnmappedChr = "*";

/**
 * Call the consensus sequence from the consensus matrix using the majority voting method.
 *
 * For each position, the function counts the frequency of each base (A, C, G, T, N, Partial, Gap) and determines
 * the majority base. If the majority base is a gap, it only calls it a gap if its frequency is strictly greater than
 * the consensus gap threshold and the forward and reverse strands agree on the base.
 *
 * If there is a tie between two bases, the lexicographic order is used to break the tie (e.g., A < C < G < T < Gap).
 *
 * @param consensus_matrix The consensus matrix containing the MSA matrix and the strand information.
 * @param options The ConsensusOptions containing the consensus parameters.
 *
 * @return A ConsensusResult struct containing the consensus sequence, quality scores, and depths.
 */
ConsensusResult MajorityVotingConsensus(const ConsensusMatrix& consensus_matrix, const ConsensusOptions& options);

// Update the consensus cluster metrics based on the clusters provided.
void UpdateConsensusMetrics(const Clusters& clusters, Metrics& metrics);

/**
 * Perform clustering and consensus on @p super_region using @p context, reading alignments and writing the consensus
 * sequences to the context's FASTQ writer. @p region_lookup_table is used to attribute each alignment to the region
 * closest to its start position so it is processed exactly once.
 */
void ClusterAndConsensusSuperRegion(ThreadContext& context,
                                    const ReadCollapserOptions& options,
                                    const SuperRegion& super_region,
                                    u32 super_region_id,
                                    const RegionLookupTable& region_lookup_table);

/**
 * Perform optimized clustering and consensus on a BAM file. This mode will:
 *  1. Perform positional clustering followed by fixed UMI clustering (if --umi-type=fixed).
 *  2. Perform provided consensus method.
 *  3. Produce consensus reads in compressed FASTQ format.
 */
void FastClusterAndConsensus(const ReadCollapserOptions& options);

}  // namespace xoos::read_collapser
