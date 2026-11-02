#include "consensus/consensus.h"

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <unordered_map>

#include <htslib/sam.h>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/alignment-reader.h>
#include <xoos/io/bam/read-group.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/read-name/consensus-read-name.h>
#include <xoos/types/int.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "clustering/cluster-downsample.h"
#include "clustering/clustering.h"
#include "consensus/base-encoder.h"
#include "consensus/consensus-matrix-builder.h"
#include "consensus/consensus-result.h"
#include "consensus/majority-voting.h"
#include "core/region-lookup.h"
#include "core/region-processor.h"
#include "core/region.h"
#include "core/thread-context.h"
#include "io/alignment-io.h"
#include "io/alignment.h"
#include "io/fastq-writer.h"
#include "io/output-path.h"
#include "metrics/metrics.h"
#include "util/duplex-util.h"
#include "util/read-util.h"
#include "util/softclip-util.h"

namespace xoos::read_collapser {

ConsensusResult MajorityVotingConsensus(const ConsensusMatrix& consensus_matrix, const ConsensusOptions& options) {
  // Initialize the consensus result with the length of the consensus sequence
  ConsensusResult consensus_result(consensus_matrix.GetConsensusLength());
  // Initialize the per-base depth and error fields in the consensus result if the aux tags are enabled. These feed
  // both the per-base array tags (ad/bd/ae/be) and the per-strand summary tags derived from them.
  if (options.include_consensus_aux_tags) {
    consensus_result.forward_per_base_depth = vec<u32>();
    consensus_result.reverse_per_base_depth = vec<u32>();
    consensus_result.forward_per_base_errors = vec<u32>();
    consensus_result.reverse_per_base_errors = vec<u32>();
  }

  // Create a majority voter to be reused for each base position (column) in the consensus matrix
  ColumnMajorityVotingWorker majority_voter(options);

  for (size_t col = 0; col < consensus_matrix.GetConsensusLength(); ++col) {
    // Compute the majority base at the current position
    ColumnConsensusResult column_result{};
    majority_voter.ComputeColumnConsensusResult(consensus_matrix, col, column_result);
    // Skip position if no reads covered it
    if (column_result.read_support == 0) {
      continue;
    }
    // Skip position if the majority call was a gap (or 'P' or 'N')
    if (column_result.majority_base == kBaseP || column_result.majority_base == kBaseN ||
        column_result.majority_base == kBaseGap) {
      continue;
    }
    // Append the final results
    consensus_result.sequence += column_result.majority_base;
    consensus_result.depths.push_back(column_result.read_support);
    consensus_result.quality_scores.push_back(column_result.qscore);

    // Record per-strand depth (bases participating in majority voting) and per-strand errors (bases disagreeing with
    // the majority base = depth - majority_count) at this position. These feed the duplex aux tags in the FASTQ
    // output when enabled. The vectors are already initialized above, so we only append here.
    if (options.include_consensus_aux_tags) {
      const auto& metrics = column_result.metrics;
      consensus_result.forward_per_base_depth->emplace_back(metrics.fwd_read_support);
      consensus_result.reverse_per_base_depth->emplace_back(metrics.rev_read_support);
      consensus_result.forward_per_base_errors->emplace_back(metrics.fwd_read_support - metrics.fwd_majority_count);
      consensus_result.reverse_per_base_errors->emplace_back(metrics.rev_read_support - metrics.rev_majority_count);
    }
  }
  return consensus_result;
}

// Update the consensus cluster metrics based on the clusters provided.
void UpdateConsensusMetrics(const Clusters& clusters, Metrics& metrics) {
  for (const auto& cluster : clusters | std::views::values) {
    const u64 cluster_size = cluster->alignments.size();
    metrics.consensus_input_reads += cluster_size;
    ++metrics.consensus_clusters;
    metrics.consensus_cluster_sizes.AddCountToHistogram(cluster_size, static_cast<u64>(1));

    bool cluster_is_full = false;
    bool cluster_is_partial = false;
    bool cluster_is_forward = false;
    bool cluster_is_reverse = false;
    for (const auto& alignment : cluster->alignments) {
      if (alignment->IsPartial()) {
        cluster_is_partial = true;
        ++metrics.consensus_input_partial_reads;
      } else {
        cluster_is_full = true;
        ++metrics.consensus_input_full_reads;
      }
      if (alignment->IsForward()) {
        cluster_is_forward = true;
      } else {
        cluster_is_reverse = true;
      }
    }
    if (cluster_is_full && cluster_is_partial) {
      ++metrics.consensus_full_and_partial_read_clusters;
    } else if (cluster_is_full) {
      ++metrics.consensus_full_read_clusters;
    } else if (cluster_is_partial) {
      ++metrics.consensus_partial_read_clusters;
    }
    if (cluster_is_forward && cluster_is_reverse) {
      ++metrics.consensus_mixed_strand_clusters;
    } else if (cluster_is_forward) {
      ++metrics.consensus_forward_strand_clusters;
    } else if (cluster_is_reverse) {
      ++metrics.consensus_reverse_strand_clusters;
    }
  }
}

namespace {

// Percentile corresponding to the median.
constexpr u8 kMedianPercentile = 50;

// Helper function to write unmapped alignments to BAM and FASTQ files.
void WriteUnmappedAlignment(const ReadCollapserOptions& options,
                            const AlignmentWriter& alignment_writer,
                            const GzipFilePtr& fastq_writer,
                            const ClusterId& cluster_id,
                            const bam1_t* const record,
                            std::string& buffer) {
  if (options.consensus_options.output_cluster_bam) {
    WriteAlignment(record, false, cluster_id, 1, alignment_writer.bam);
  }
  u32 avg_cluster_size = 1;
  if (options.consensus_options.duplex_library_type != HDDeconvolutionType::kNone) {
    UpdateUnmappedDeconvolvedQualities(record);
    avg_cluster_size = GetUnmappedDeconvolvedAvgClusterSize(record);
  }
  const auto seq = GetSequence(bam_get_seq(record), 0, record->core.l_qseq);
  const auto qualities = GetQualities(bam_get_qual(record), 0, record->core.l_qseq);
  // Cluster ID followed by the number of partial fwd, partial rev, full fwd, and full rev reads (which will
  // always be 0) and the average cluster size
  const std::string read_name = read_name::FormatConsensusReadName(
      cluster_id.super_region_id, cluster_id.count, {.avg_family_size = avg_cluster_size});
  GzipWriteFastq(fastq_writer.get(), seq, qualities, read_name, buffer);
}

void CallConsensusOnLeftSoftclipsAndPrependResult(const vec<AlignmentPtr>& reads_in_cluster,
                                                  const ReadCollapserOptions& options,
                                                  ConsensusResult& consensus_result) {
  const vec<SoftclipSequence> left_softclips = GetLeftSoftclipSequences(reads_in_cluster);
  const bool has_left_softclips =
      std::ranges::any_of(left_softclips, [](const SoftclipSequence& seq) { return !seq.sequence.empty(); });
  if (has_left_softclips) {
    const ConsensusMatrix left_softclip_matrix = BuildSoftclipConsensusMatrix(
        left_softclips, reads_in_cluster, options.consensus_options.duplex_library_type != HDDeconvolutionType::kNone);
    consensus_result = MajorityVotingConsensus(left_softclip_matrix, options.consensus_options) + consensus_result;
  }
}

void CallConsensusOnRightSoftclipsAndAppendResult(const vec<AlignmentPtr>& reads_in_cluster,
                                                  const ReadCollapserOptions& options,
                                                  ConsensusResult& consensus_result) {
  const vec<SoftclipSequence> right_softclips = GetRightSoftclipSequences(reads_in_cluster);
  const bool has_right_softclips =
      std::ranges::any_of(right_softclips, [](const SoftclipSequence& seq) { return !seq.sequence.empty(); });
  if (has_right_softclips) {
    const ConsensusMatrix right_softclip_matrix = BuildSoftclipConsensusMatrix(
        right_softclips, reads_in_cluster, options.consensus_options.duplex_library_type != HDDeconvolutionType::kNone);
    consensus_result = consensus_result + MajorityVotingConsensus(right_softclip_matrix, options.consensus_options);
  }
}

/**
 * Call consensus on the clusters, update consensus metrics, and write the results to the output FASTQ file.
 *
 * This function handles unmapped regions by writing the reads directly to the FASTQ file.
 * For mapped regions, it performs clustering, deconvolution (if enabled), and consensus calling
 * before writing the results to the FASTQ file.
 *
 * @param options Read collapser options.
 * @param clusters Clusters containing the alignments to process.
 * @param fastq_writer Gzip file pointer for writing the output FASTQ file.
 * @param metrics The metrics to accumulate consensus metrics into.
 */
// NOLINTNEXTLINE(readability-function-cognitive-complexity) - sequential consensus calling and result writing
void CallConsensusAndWriteResults(const ReadCollapserOptions& options,
                                  const Clusters& clusters,
                                  const GzipFilePtr& fastq_writer,
                                  Metrics& metrics) {
  Clusters consensus_clusters;
  for (const auto& [cluster_id, cluster] : clusters) {
    // Filter out small clusters
    if (cluster->alignments.size() < options.consensus_options.min_cluster_size) {
      ++metrics.consensus_discarded_by_size_clusters;
      continue;
    }
    consensus_clusters[cluster_id] = cluster;
  }
  // Subsample large clusters to the maximum number of reads per cluster before collecting metrics
  // and calling consensus
  for (const auto& [cluster_id, cluster] : consensus_clusters) {
    if (cluster->alignments.size() > options.consensus_options.max_cluster_size) {
      const auto num_reads_before_downsampling = cluster->alignments.size();
      DownsampleReadsInCluster(cluster->alignments, options.consensus_options.max_cluster_size, cluster_id);
      const auto num_reads_after_downsampling = cluster->alignments.size();
      metrics.consensus_discarded_by_subsampling_reads +=
          (num_reads_before_downsampling - num_reads_after_downsampling);
    }
  }
  UpdateConsensusMetrics(consensus_clusters, metrics);
  // Reused across clusters and calls so each GzipWriteFastq call does not reallocate its output buffer.
  thread_local std::string fastq_buffer;
  for (const auto& [cluster_id, cluster] : consensus_clusters) {
    const bool hd_deconvolve_enabled = (options.consensus_options.duplex_library_type != HDDeconvolutionType::kNone);
    bool handle_softclips = !options.consensus_options.skip_softclips;
    // If trimming overhangs (soft-clip handling) would produce an empty reference span
    // (i.e. left boundary >= right boundary), disable soft-clip handling and keep overhangs
    // to avoid generating empty consensus sequences.
    if (handle_softclips) {
      const u32 leftmost_rpos = LeftmostRefPos(cluster->alignments, true);
      const u32 rightmost_rpos = RightmostRefPos(cluster->alignments, true);
      if (leftmost_rpos >= rightmost_rpos) {
        handle_softclips = false;
      }
    }
    const ConsensusMatrix consensus_matrix =
        BuildConsensusMatrix(cluster->alignments, hd_deconvolve_enabled, handle_softclips);
    ConsensusResult consensus_result = MajorityVotingConsensus(consensus_matrix, options.consensus_options);
    // If softclip handling is enabled, we also need to create two additional consensus matrices
    // one for the left softclip and one for the right softclip
    // and merge the softclip consensus results with the main consensus result
    if (handle_softclips) {
      CallConsensusOnLeftSoftclipsAndPrependResult(cluster->alignments, options, consensus_result);
      CallConsensusOnRightSoftclipsAndAppendResult(cluster->alignments, options, consensus_result);
    }

    // Trim the ends of the consensus result to remove low coverage bases
    consensus_result.TrimEnds(options.consensus_options.min_trim_read_support);

    if (options.consensus_options.min_consensus_read_length.has_value() &&
        consensus_result.sequence.length() < options.consensus_options.min_consensus_read_length.value()) {
      // Discard consensus results that are too short
      ++metrics.consensus_discarded_by_length_reads;
      continue;
    }

    // Write the consensus result to the output FASTQ file
    if (!consensus_result.sequence.empty()) {
      const ClusterMetadata metadata = consensus_matrix.GetMetadata();
      const std::string read_name =
          read_name::FormatConsensusReadName(cluster_id.super_region_id,
                                             cluster_id.count,
                                             {.forward_partial = metadata.forward_partial,
                                              .reverse_partial = metadata.reverse_partial,
                                              .forward_full = metadata.forward_full,
                                              .reverse_full = metadata.reverse_full,
                                              .avg_family_size = consensus_result.MeanClusterSize()});
      const bool has_forward = (metadata.forward_partial + metadata.forward_full) > 0;
      const bool has_reverse = (metadata.reverse_partial + metadata.reverse_full) > 0;
      const bool is_same_strand = !has_forward || !has_reverse;
      const bool meets_same_strand_depth_threshold =
          is_same_strand &&
          consensus_result.MeanClusterSize() >= options.consensus_options.min_same_strand_cluster_size;
      const bool meets_mixed_strand_depth_threshold =
          !is_same_strand &&
          consensus_result.MeanClusterSize() >= options.consensus_options.min_mixed_strand_cluster_size;
      if (meets_same_strand_depth_threshold || meets_mixed_strand_depth_threshold) {
        ++metrics.total_consensus_reads;
        // Only emit per-base arrays for strands present in the cluster: a same-strand cluster's opposite-strand
        // ad/ae or bd/be arrays would be uniformly zero and carry no information.
        GzipWriteFastq(fastq_writer.get(), consensus_result, read_name, fastq_buffer, has_forward, has_reverse);
      } else if (is_same_strand) {
        ++metrics.consensus_discarded_by_same_strand_cluster_size_reads;
      } else {
        ++metrics.consensus_discarded_by_mixed_strand_cluster_size_reads;
      }
    }
  }
}

}  // namespace

// Perform clustering and consensus on a super region
void ClusterAndConsensusSuperRegion(ThreadContext& context,
                                    const ReadCollapserOptions& options,
                                    const SuperRegion& super_region,
                                    const u32 super_region_id,
                                    const RegionLookupTable& region_lookup_table) {
  const auto& fastq_writer = context.FastqWriter();
  auto& metrics = context.GetMetrics();
  u32 clusters_id = 0;
  auto create_cluster_id = [&super_region_id, &clusters_id]() -> ClusterId {
    return {super_region_id, clusters_id++};  // NOLINT
  };
  // Per-batch action: cluster the alignments, optionally emit a cluster BAM, then call consensus and write the results.
  const auto process_batch =
      [&options, &create_cluster_id, &metrics, &context, &fastq_writer](vec<AlignmentPtr>& alignments) {
        const Clusters clusters = ClusterAlignments(options, alignments, create_cluster_id, metrics);
        if (options.consensus_options.output_cluster_bam) {
          // Emit clustered alignments to the cluster BAM, always including cluster info.
          WriteAlignments(alignments, false, true, context.BamWriter().bam);
        }
        CallConsensusAndWriteResults(options, clusters, fastq_writer, metrics);
      };
  for (size_t sub_region_idx = 0; sub_region_idx < super_region.size(); ++sub_region_idx) {
    const auto& region = super_region.at(sub_region_idx);
    // Ownership strategy: closest-start. An alignment is owned by the region whose start is closest to the
    // alignment's start, preserving cluster integrity at region boundaries.
    const auto owns_read =
        [&region_lookup_table, &region, &super_region_id, &sub_region_idx](const bam1_t* const record) {
          const auto closest_region = region_lookup_table.FindClosestOverlappingRegion(
              region.tid,
              {.start = record->core.pos, .end = bam_endpos(record)},
              {.super_region_id = super_region_id, .sub_region_index = sub_region_idx});
          return closest_region.has_value() && std::cmp_equal(closest_region->start, region.start) &&
                 std::cmp_equal(closest_region->end, region.end);
        };
    ProcessRegionInBatches(context, options, region, owns_read, process_batch);
  }
}

void FastClusterAndConsensus(const ReadCollapserOptions& options) {
  ValidateOutputFilesDoNotExist(options, ReadCollapserMode::kConsensus);

  // For consensus, we do not need to limit the number of super regions because the output does not need to be sorted.
  auto super_regions = DetermineSuperRegions(options);

  auto super_region_count = super_regions.size();
  Logging::Info("Processing {} super regions", super_region_count);

  // One context per thread; each owns its own reader, FASTQ writer, optional cluster BAM writer, and metrics.
  // Do not write an index for the cluster BAM because, unlike mark duplicate, reads may be processed out of order for
  // consensus.
  std::vector<ThreadContext> contexts;
  contexts.reserve(options.threads);
  for (size_t i = 0; i < options.threads; ++i) {
    auto& context = contexts.emplace_back(options);
    if (options.consensus_options.output_cluster_bam) {
      context.OpenBamWriter(ConsensusClusterBamPath(options, i), /*write_index=*/false);
    }
    context.OpenFastqWriter(ConsensusFastqPath(options, i), options.consensus_options.compression_level);
  }

  // Create a region lookup table for finding the closest region to the start position of an alignment
  const auto region_lookup_table = RegionLookupTable(super_regions);

  tf::Taskflow taskflow;
  tf::Executor executor{options.threads};
  for (u32 super_region_id = 0; super_region_id < super_region_count; ++super_region_id) {
    // Need to capture super_region_id by value because by the time the task is
    // actually executed, it would be out of scope
    auto& super_region = super_regions.at(super_region_id);
    taskflow.emplace([&options, super_region_id, &super_region, &contexts, &executor, &region_lookup_table] {
      ThreadContext& context = contexts.at(static_cast<size_t>(executor.this_worker_id()));
      try {
        // For consensus, each super region only contains subregions on the same tid
        // so we can check the first subregion to determine if it is the special unmapped region
        if (!super_region.empty() && super_region.front().tid == HTS_IDX_NOCOOR) {
          u32 clusters_id = 0;
          auto create_cluster_id = [&super_region_id, &clusters_id]() -> ClusterId {
            return {super_region_id, clusters_id++};  // NOLINT
          };
          // This is the super region for unmapped reads. Write unmapped reads directly to output files.
          // Reused across records so each GzipWriteFastq call does not reallocate its output buffer.
          // thread_local (static storage) so the lambda accesses it directly rather than capturing it.
          static thread_local std::string fastq_buffer;
          auto writer = [&options, &context](const ClusterId& cluster_id, const bam1_t* const record) {
            WriteUnmappedAlignment(
                options, context.BamWriter(), context.FastqWriter(), cluster_id, record, fastq_buffer);
          };
          ReadAndWriteUnmappedReads(options, create_cluster_id, context.Reader(), writer, context.GetMetrics());
        } else {
          ClusterAndConsensusSuperRegion(context, options, super_region, super_region_id, region_lookup_table);
        }
      } catch (const std::exception& e) {
        if (!super_region.empty() && super_region.front().tid == HTS_IDX_NOCOOR) {
          Logging::Error("Error processing unmapped reads: {}", e.what());
        } else {
          Logging::Error("Error processing super region {}:{}-{}: {}",
                         io::SamHdrTid2Name(context.Reader().hdr.get(), super_region.front().tid),
                         super_region.front().start,
                         super_region.back().end,
                         e.what());
        }
        throw;
      }
    });
  }
  executor.run(taskflow).get();

  // Merge per-context metrics into a single set.
  vec<Metrics> worker_metrics;
  worker_metrics.reserve(contexts.size());
  for (const auto& context : contexts) {
    worker_metrics.emplace_back(context.GetMetrics());
  }
  auto metrics_total = SumMetrics(worker_metrics);
  // update the consensus median cluster size now that all data is available
  if (!metrics_total.consensus_cluster_sizes.IsEmpty()) {
    metrics_total.consensus_median_cluster_size =
        histogram::ComputePercentile(metrics_total.consensus_cluster_sizes, kMedianPercentile);
  }
  const io::TsvMetadata tsv_metadata{.command_line_info = options.command_line_info,
                                     .read_group = io::GetFirstReadGroupFromBam(options.bam_input)};
  metrics_total.WriteAllMetricsToTsv(ReadCollapserMode::kConsensus, options, tsv_metadata);

  Logging::Info("Done");
}

}  // namespace xoos::read_collapser
