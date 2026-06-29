#pragma once

#include <xoos/types/float.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "consensus/consensus-matrix.h"
#include "io/alignment.h"

namespace xoos::read_collapser {

using InsertionMap = std::unordered_map<s64, vec<std::string>>;

/**
 * Build a consensus matrix for the cluster that can be used for calling the consensus sequence.
 * A multiple sequence alignment of the sequences is produced by aligning the insertions in the sequences
 * and using the insertions as anchors to construct the MSA.
 *
 * @param reads_in_cluster A vector of AlignmentPtr representing the reads in the cluster.
 * @param hd_deconvolve_enabled Whether to populate homopolymer information in the consensus matrix. Homopolymer
 * information is only necessary if HD deconvolution is enabled.
 * @param trim_overhangs Whether to trim overhangs when building the consensus matrix. If this is set to true,
 * reads with soft-clips will have their overhang bases trimmed so that their start and end positions are
 * consistent (i.e. all reads with soft clipped bases will have "blunt" ends). The trimmed bases can then
 * be handled separately along with the soft-clipped bases.
 *
 * @return A struct containing a 2D vector representing the consensus matrix as well as another vector
 * indicating whether each read is on the forward or reverse strand. The consensus matrix is a 2D vector
 * where each row corresponds to a read in the cluster and each column corresponds to a position in the
 * consensus sequence. The matrix is filled with the bases from the original sequences, the aligned
 * insertions, and gaps ('-') where appropriate.
 */
ConsensusMatrix BuildConsensusMatrix(const std::vector<AlignmentPtr>& reads_in_cluster,
                                     bool hd_deconvolve_enabled,
                                     bool trim_overhangs);

/**
 * Given a map of aligned insertions and the cluster from which the insertions are found,
 * build a consensus matrix for the cluster that can be used for calling the consensus sequence.
 *
 * @param aligned_insertions A map where the key is the reference position and the value is a vector of strings
 * representing the gapped aligned insertions at that position for each read in the cluster.
 * @param reads_in_cluster A vector of AlignmentPtr representing the reads in the cluster.
 * @param hd_deconvolve_enabled Whether to populate homopolymer information in the consensus matrix. Homopolymer
 * information is only necessary if HD deconvolution is enabled.
 * @param trim_overhangs Whether to trim overhangs when building the consensus matrix. If this is set to true,
 * reads with soft-clips will have their overhang bases trimmed so that their start and end positions are
 * consistent (i.e. all reads with soft clipped bases will have "blunt" ends). The trimmed bases can then
 * be handled separately along with the soft-clipped bases.
 *
 * @return A struct containing a 2D vector representing the consensus matrix as well as another vector
 * indicating whether each read is on the forward or reverse strand. The consensus matrix is a 2D vector
 * where each row corresponds to a read in the cluster and each column corresponds to a position in the
 * consensus sequence. The matrix is filled with the bases from the original sequences, the aligned
 * insertions, and gaps ('-') where appropriate.
 */
ConsensusMatrix BuildConsensusMatrix(const InsertionMap& aligned_insertions,
                                     const vec<AlignmentPtr>& reads_in_cluster,
                                     bool hd_deconvolve_enabled,
                                     bool trim_overhangs);

/**
 * Build a consensus matrix from a set of sequences without alignment information. Since we cannot use insertions
 * as anchors to build the MSA, we perform a full MSA on the sequences to when building the consensus matrix.
 *
 * If the input sequences are empty, return an empty consensus matrix.
 *
 * @param sequences A vector of strings representing the sequences in the cluster.
 * @param reads_in_cluster A vector of AlignmentPtr representing the reads in the cluster. Note that we do not
 * use the sequences from the alignments, but we do extract strand and duplex strand information from the alignments.
 */
ConsensusMatrix BuildConsensusMatrix(const vec<std::string>& sequences, const vec<AlignmentPtr>& reads_in_cluster);

/**
 * @brief Given a vector of alignments in a cluster, find the insertions each alignment and align them.
 * @param reads_in_cluster A vector of AlignmentPtr representing the reads in the cluster.
 * @return A map where the key is the reference position where at least one read in the cluster has an insertion,
 * and the value is a vector of strings representing the gapped aligned insertions at that position for each read in the
 * cluster. Leading insertions (i.e., insertions at the beginning of a read or immediately following a left soft clip)
 * and trailing insertions (i.e. insertions at the end of a read or immediately preceding a right soft clip) are
 * ignored.
 */
InsertionMap GetAlignedInsertions(const vec<AlignmentPtr>& reads_in_cluster);

}  // namespace xoos::read_collapser
