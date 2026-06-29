#include "consensus/consensus.h"

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <string_view>
#include <unordered_map>

#include <htslib/sam.h>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "clustering/cluster-downsample.h"
#include "clustering/clustering.h"
#include "consensus/base-encoder.h"
#include "consensus/consensus-matrix-builder.h"
#include "consensus/consensus-result.h"
#include "consensus/majority-voting.h"
#include "core/region-lookup.h"
#include "core/region.h"
#include "io/alignment-io.h"
#include "io/alignment.h"
#include "io/fastq-writer.h"
#include "metrics/metrics.h"
#include "util/duplex-util.h"
#include "util/read-util.h"
#include "util/softclip-util.h"

namespace xoos::read_collapser {

ConsensusResult MajorityVotingConsensus(const ConsensusMatrix& consensus_matrix, const ReadCollapserOptions& options) {
  // Initialize the consensus result with the length of the consensus sequence
  ConsensusResult consensus_result(consensus_matrix.GetConsensusLength());
  // Initialize the per-base depth and majority count fields in the consensus result if the options are enabled
  if (options.include_per_base_read_support_tags) {
    consensus_result.forward_per_base_depth = vec<u32>();
    consensus_result.reverse_per_base_depth = vec<u32>();
  }
  if (options.include_per_base_majority_count_tags) {
    consensus_result.forward_per_base_majority_count = vec<u32>();
    consensus_result.reverse_per_base_majority_count = vec<u32>();
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

    // Record the base counts from the forward and reverse strands, as well as the number of forward and reverse bases
    // that agree with the majority base at this position. These will be used to encode the per-base depth tags in the
    // FASTQ output if the option is enabled.
    if (options.include_per_base_read_support_tags) {
      // If the option is enabled, then the vectors should already be initialized
      // so we just need to append the values
      consensus_result.forward_per_base_depth->emplace_back(column_result.metrics.fwd_read_support);
      consensus_result.reverse_per_base_depth->emplace_back(column_result.metrics.rev_read_support);
    }
    if (options.include_per_base_majority_count_tags) {
      // If the option is enabled, then the vectors should already be initialized
      // so we just need to append the values
      consensus_result.forward_per_base_majority_count->emplace_back(column_result.metrics.fwd_majority_count);
      consensus_result.reverse_per_base_majority_count->emplace_back(column_result.metrics.rev_majority_count);
    }
  }
  return consensus_result;
}

// Update the consensus cluster metrics based on the clusters provided.
void UpdateConsensusMetrics(const Clusters& clusters) {
  auto& metrics = ConcurrentMetrics::Get();
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

// Helper function to write unmapped alignments to BAM and FASTQ files.
static void WriteUnmappedAlignment(const ReadCollapserOptions& options,
                                   const AlignmentWriter& alignment_writer,
                                   const GzipFilePtr& fastq_writer,
                                   const ClusterId& cluster_id,
                                   const bam1_t* record) {
  if (options.output_cluster_bam) {
    WriteAlignment(record, false, cluster_id, 1, alignment_writer.bam);
  }
  u32 avg_cluster_size = 1;
  if (options.duplex_library_type != HDDeconvolutionType::kNone) {
    UpdateUnmappedDeconvolvedQualities(record);
    avg_cluster_size = GetUnmappedDeconvolvedAvgClusterSize(record);
  }
  const auto seq = GetSequence(bam_get_seq(record), 0, record->core.l_qseq);
  const auto qualities = GetQualities(bam_get_qual(record), 0, record->core.l_qseq);
  // Cluster ID followed by the number of partial fwd, partial rev, full fwd, and full rev reads (which will
  // always be 0) and the average cluster size
  const std::string read_name =
      fmt::format("{}:{}-0-0-0-0-{}", cluster_id.super_region_id, cluster_id.count, avg_cluster_size);
  GzipWriteFastq(fastq_writer.get(), seq, qualities, read_name);
}

static void CallConsensusOnLeftSoftclipsAndPrependResult(const vec<AlignmentPtr>& reads_in_cluster,
                                                         const ReadCollapserOptions& options,
                                                         ConsensusResult& consensus_result) {
  const vec<std::string> left_softclips = GetLeftSoftclipSequences(reads_in_cluster);
  const bool has_left_softclips =
      std::ranges::any_of(left_softclips, [](const std::string_view seq) { return !seq.empty(); });
  if (has_left_softclips) {
    const ConsensusMatrix left_softclip_matrix = BuildConsensusMatrix(left_softclips, reads_in_cluster);
    consensus_result = MajorityVotingConsensus(left_softclip_matrix, options) + consensus_result;
  }
}

static void CallConsensusOnRightSoftclipsAndAppendResult(const vec<AlignmentPtr>& reads_in_cluster,
                                                         const ReadCollapserOptions& options,
                                                         ConsensusResult& consensus_result) {
  const vec<std::string> right_softclips = GetRightSoftclipSequences(reads_in_cluster);
  const bool has_right_softclips =
      std::ranges::any_of(right_softclips, [](const std::string_view seq) { return !seq.empty(); });
  if (has_right_softclips) {
    const ConsensusMatrix right_softclip_matrix = BuildConsensusMatrix(right_softclips, reads_in_cluster);
    consensus_result = consensus_result + MajorityVotingConsensus(right_softclip_matrix, options);
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
 */
static void CallConsensusAndWriteResults(const ReadCollapserOptions& options,
                                         const Clusters& clusters,
                                         const GzipFilePtr& fastq_writer) {
  auto& metrics = ConcurrentMetrics::Get();
  Clusters consensus_clusters;
  for (const auto& [cluster_id, cluster] : clusters) {
    // Filter out small clusters
    if (cluster->alignments.size() < options.min_cluster_size) {
      ++metrics.consensus_discarded_by_size_clusters;
      continue;
    }
    consensus_clusters[cluster_id] = cluster;
  }
  // Subsample large clusters to the maximum number of reads per cluster before collecting metrics
  // and calling consensus
  for (const auto& [cluster_id, cluster] : consensus_clusters) {
    if (cluster->alignments.size() > options.max_cluster_size) {
      const auto num_reads_before_downsampling = cluster->alignments.size();
      DownsampleReadsInCluster(cluster->alignments, options.max_cluster_size);
      const auto num_reads_after_downsampling = cluster->alignments.size();
      metrics.consensus_discarded_by_subsampling_reads +=
          (num_reads_before_downsampling - num_reads_after_downsampling);
    }
  }
  UpdateConsensusMetrics(consensus_clusters);
  for (const auto& [cluster_id, cluster] : consensus_clusters) {
    const bool hd_deconvolve_enabled = (options.duplex_library_type != HDDeconvolutionType::kNone);
    const bool handle_softclips = options.include_softclips;
    const ConsensusMatrix consensus_matrix =
        BuildConsensusMatrix(cluster->alignments, hd_deconvolve_enabled, handle_softclips);
    ConsensusResult consensus_result = MajorityVotingConsensus(consensus_matrix, options);
    // If softclip handling is enabled, we also need to create two additional consensus matrices
    // one for the left softclip and one for the right softclip
    // and merge the softclip consensus results with the main consensus result
    if (handle_softclips) {
      CallConsensusOnLeftSoftclipsAndPrependResult(cluster->alignments, options, consensus_result);
      CallConsensusOnRightSoftclipsAndAppendResult(cluster->alignments, options, consensus_result);
    }

    // Trim the ends of the consensus result to remove low coverage bases
    consensus_result.TrimEnds(options.min_trim_read_support);

    if (options.min_consensus_read_length.has_value() &&
        consensus_result.sequence.length() < options.min_consensus_read_length.value()) {
      // Discard consensus results that are too short
      ++metrics.consensus_discarded_by_length_reads;
      continue;
    }

    // Write the consensus result to the output FASTQ file
    if (!consensus_result.sequence.empty()) {
      const ClusterMetadata metadata = consensus_matrix.GetMetadata();
      const std::string read_name = fmt::format("{}:{}-{}-{}-{}-{}-{}",
                                                cluster_id.super_region_id,
                                                cluster_id.count,
                                                metadata.forward_partial,
                                                metadata.reverse_partial,
                                                metadata.forward_full,
                                                metadata.reverse_full,
                                                consensus_result.MeanClusterSize());
      const bool has_forward = (metadata.forward_partial + metadata.forward_full) > 0;
      const bool has_reverse = (metadata.reverse_partial + metadata.reverse_full) > 0;
      const bool is_same_strand = !has_forward || !has_reverse;
      const bool meets_same_strand_depth_threshold =
          is_same_strand && consensus_result.MeanClusterSize() >= options.min_same_strand_cluster_size;
      const bool meets_mixed_strand_depth_threshold =
          !is_same_strand && consensus_result.MeanClusterSize() >= options.min_mixed_strand_cluster_size;
      if (meets_same_strand_depth_threshold || meets_mixed_strand_depth_threshold) {
        ++metrics.total_consensus_reads;
        GzipWriteFastq(fastq_writer.get(), consensus_result, read_name);
      } else if (is_same_strand) {
        ++metrics.consensus_discarded_by_same_strand_cluster_size_reads;
      } else {
        ++metrics.consensus_discarded_by_mixed_strand_cluster_size_reads;
      }
    }
  }
}

// Perform clustering and consensus on a super region
void ClusterAndConsensusSuperRegion(const ReadCollapserOptions& options,
                                    const u32 super_region_id,
                                    const SuperRegion& super_region,
                                    const AlignmentReader& alignment_reader,
                                    const AlignmentWriter& alignment_writer,
                                    const GzipFilePtr& fastq_writer,
                                    const RegionLookupTable& region_lookup_table) {
  auto& metrics = ConcurrentMetrics::Get();
  u32 clusters_id = 0;
  auto create_cluster_id = [&super_region_id, &clusters_id]() -> ClusterId {
    return {super_region_id, clusters_id++};  // NOLINT
  };
  for (size_t sub_region_idx = 0; sub_region_idx < super_region.size(); ++sub_region_idx) {
    const auto& region = super_region.at(sub_region_idx);
    const auto itr = io::SamItrQueryI(alignment_reader.idx.get(), region.tid, region.start, region.end);
    vec<AlignmentPtr> alignments;
    // Read alignments for the region in batches
    // If the batch size is reached, cluster the alignments and call consensus before continuing to the next batch
    while (true) {
      auto record = io::Bam1Ptr{bam_init1()};
      const bool has_more_reads = io::SamItrNext(alignment_reader.bam.get(), itr.get(), record.get());
      const bool reached_batch_size_limit = (alignments.size() >= options.batch_size);
      if ((!has_more_reads || reached_batch_size_limit) && !alignments.empty()) {
        Clusters clusters;
        // Do clustering on the alignments and update metrics
        clusters = ClusterAlignments(options, alignments, create_cluster_id);
        if (options.output_cluster_bam) {
          // Output clustered alignments to BAM file if specified
          // Always include cluster info in the BAM output
          WriteAlignments(alignments, false, true, alignment_writer.bam);
        }
        CallConsensusAndWriteResults(options, clusters, fastq_writer);
        alignments.clear();  // Clear the alignments for the next batch
      }
      if (!has_more_reads) {
        break;  // No more records to read
      }
      // Check if the start position of the current region is the closest to the start position of the alignment
      // If not, skip the alignment as an alignment should only be handled by the region closest to its start position
      // to avoid processing the same alignment multiple times.
      const auto closest_region = region_lookup_table.FindClosestOverlappingRegion(
          region.tid, {record->core.pos, bam_endpos(record.get())}, {super_region_id, sub_region_idx});
      const bool current_region_is_closest =
          (closest_region.has_value() && std::cmp_equal(closest_region->start, region.start) &&
           std::cmp_equal(closest_region->end, region.end));
      if (!current_region_is_closest) {
        // If the current region is not the closest to the start position of the alignment, skip the alignment
        continue;
      }
      if (ShouldDiscardAlignment(options, record.get(), metrics)) {
        continue;
      }
      alignments.emplace_back(std::make_shared<Alignment>(std::move(record)));
    }
  }
}

void FastClusterAndConsensus(const ReadCollapserOptions& options) {
  // Reset Metrics instance so that test cases do not share the same data (due to the static Metrics instance member)
  ConcurrentMetrics::Reset();

  // For consensus, we do not need to limit the number of super regions because the output does not need to be sorted.
  auto super_regions = DetermineSuperRegions(options);

  auto super_region_count = super_regions.size();
  Logging::Info("Processing {} super regions", super_region_count);

  auto alignment_readers = OpenAlignmentReaders(options.bam_input, options.threads);

  // Only output cluster BAM if options.output_cluster_bam is true
  // Otherwise, create empty AlignmentWriter objects to avoid bad access during taskflow processing
  auto alignment_writers = vec<AlignmentWriter>();
  for (size_t i = 0; i < options.threads; ++i) {
    if (options.output_cluster_bam) {
      const auto bam_filename = options.output_dir / fmt::format("output.{:04}.bam", i);
      // Do not write index file for BAM output because unlike mark duplicate, reads may be processed out of order for
      // consensus
      const auto hdr = alignment_readers.at(i).hdr.get();
      io::SamHdrAddPgLine(
          hdr, io::PgHdrLine{options.program_name, options.program_name, options.version, options.command_line});
      alignment_writers.emplace_back(OpenAlignmentWriter(bam_filename, hdr, false));
    } else {
      alignment_writers.emplace_back(AlignmentWriter{});
    }
  }

  auto fastq_writers = vec<GzipFilePtr>();
  for (size_t i = 0; i < options.threads; ++i) {
    // TODO: Should create an uncompressed (.gz) output file if compression level is 0
    fs::path output_fastq = options.output_dir / fmt::format("output.{:04}.fastq.gz", i);
    fastq_writers.emplace_back(OpenGzipFile(output_fastq, options.compression_level));
  }

  // Create a region lookup table for finding the closest region to the start position of an alignment
  const auto region_lookup_table = RegionLookupTable(super_regions);

  tf::Taskflow taskflow;
  tf::Executor executor{options.threads};
  for (u32 super_region_id = 0; super_region_id < super_region_count; ++super_region_id) {
    std::function<void()> super_region_func;
    // Need to capture super_region_id by value because by the time the task is
    // actually executed, it would be out of scope
    auto& super_region = super_regions.at(super_region_id);
    super_region_func = [&options,
                         super_region_id,
                         &super_region,
                         &alignment_readers,
                         &alignment_writers,
                         &fastq_writers,
                         &executor,
                         &region_lookup_table] {
      try {
        const AlignmentReader& alignment_reader = alignment_readers.at(executor.this_worker_id());
        const AlignmentWriter& alignment_writer = alignment_writers.at(executor.this_worker_id());
        const GzipFilePtr& fastq_writer = fastq_writers.at(executor.this_worker_id());
        // For consensus, each super region only contains subregions on the same tid
        // so we can check the first subregion to determine if it is the special unmapped region
        if (!super_region.empty() && super_region.front().tid == HTS_IDX_NOCOOR) {
          u32 clusters_id = 0;
          auto create_cluster_id = [&super_region_id, &clusters_id]() -> ClusterId {
            return {super_region_id, clusters_id++};  // NOLINT
          };
          // This is the super region for unmapped reads. Write unmapped reads directly to output files.
          auto writer = [&options, &alignment_writer, &fastq_writer](const ClusterId& cluster_id,
                                                                     const bam1_t* const record) {
            WriteUnmappedAlignment(options, alignment_writer, fastq_writer, cluster_id, record);
          };
          ReadAndWriteUnmappedAlignments(options, create_cluster_id, alignment_reader, writer);
        } else {
          ClusterAndConsensusSuperRegion(options,
                                         super_region_id,
                                         super_region,
                                         alignment_reader,
                                         alignment_writer,
                                         fastq_writer,
                                         region_lookup_table);
        }
      } catch (const std::exception& e) {
        if (!super_region.empty() && super_region.front().tid == HTS_IDX_NOCOOR) {
          Logging::Error("Error processing unmapped reads: {}", e.what());
        } else {
          Logging::Error("Error processing super region {}:{}-{}: {}",
                         io::SamHdrTid2Name(alignment_readers.at(0).hdr.get(), super_region.front().tid),
                         super_region.front().start,
                         super_region.back().end,
                         e.what());
        }
        throw;
      }
    };
    taskflow.emplace(super_region_func);
  }
  executor.run(taskflow).get();

  auto metrics_total = ConcurrentMetrics::GetTotal();
  // update the consensus median cluster size now that all data is available
  if (!metrics_total.consensus_cluster_sizes.IsEmpty()) {
    metrics_total.consensus_median_cluster_size =
        histogram::ComputePercentile(metrics_total.consensus_cluster_sizes, 50);
  }
  metrics_total.WriteAllMetricsToTsv(
      options.output_dir, ReadCollapserMode::kConsensus, options.cluster_by_umi, options.version, options.command_line);

  Logging::Info("Done");
}

}  // namespace xoos::read_collapser
