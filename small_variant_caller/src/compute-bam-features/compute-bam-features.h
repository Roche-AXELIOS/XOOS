#pragma once

#include <atomic>
#include <memory>

#include <xoos/io/bed-region.h>

#include "compute-bam-region-features.h"
#include "core/config.h"
#include "progress-meter.h"
#include "ref-region.h"

namespace xoos::svc {

// consolidated parameters from CLI options and config JSON
struct ComputeBamFeaturesCliParams {
  vec<fs::path> bam_input{};
  fs::path genome{};
  std::optional<fs::path> output_file{};
  SVCConfig config{};
  u8 min_mapq{};
  u8 min_bq{};
  float min_allowed_distance_from_end{};
  u32 min_family_size{};
  std::optional<u32> max_read_variant_count{};
  std::optional<float> max_read_variant_count_normalized{};
  size_t threads{};
  std::optional<vec<io::BedRegion>> bed_regions{};
  Workflow workflow{};
  bool filter_homopolymer{};
  u64 min_homopolymer_length{};
  bool duplex{};
  std::optional<std::string> tumor_read_group{};
  u64 max_region_size_per_thread{};
  std::optional<fs::path> skip_variants_vcf{};
  bool decode_yc{};
  yc_decode::BaseType min_base_type{};
};

/**
 * @brief Compute features, but on multiple threads as determined by the BED regions
 * @param cli_params A ParallelComputeFeaturesParam struct that contains user specified options
 * @returns The number of bases covered by reads processed
 */
u32 ParallelComputeBamFeatures(const ComputeBamFeaturesCliParams& cli_params);

/**
 * @brief Class responsible for orchestrating the computation of features from BAM files,
 *        including parallelization, region partitioning, and feature extraction.
 */
class ComputeBamFeatures {
 public:
  /**
   * @brief Compute features, but on multiple threads as determined by the BED regions
   * @param cli_params A ParallelComputeFeaturesParam struct that contains user specified options
   * @returns The number of bases covered by reads processed
   */
  u32 ParallelComputeBamFeatures(const ComputeBamFeaturesCliParams& cli_params);

  /**
   * @brief Given a list of BAM files, this function reads the headers and indices to compute a
   *        consolidated map of RefRegion structs that store information relating to the contigs
   *        contained within the BAM headers.
   * @param bam_paths A vector of paths to BAM files
   * @return A map of chromosome names to their corresponding RefRegion structs
   */
  StrMap<RefRegion> GetAllRefRegions(const vec<fs::path>& bam_paths);

  /**
   * @brief Loads reference sequences from the genome file for chromosomes that will be processed.
   *        Only loads sequences for chromosomes present in BED regions (if specified) to optimize memory usage.
   * @param genome Path to the reference genome FASTA file
   * @param bed_regions Map of chromosome to intervals from BED file (if provided)
   * @param contig_map Map of reference regions containing alignment coverage information
   * @return Map of chromosome names to their reference sequences
   */
  void LoadReferenceSequences(const fs::path& genome,
                              const StrMap<vec<Interval>>& bed_regions,
                              const StrMap<RefRegion>& contig_map);

  /**
   * @brief Sets up workers and alignment readers for parallel feature extraction.
   * @param cli_params A ParallelComputeFeaturesParam struct that contains user specified options
   * @param params A ComputeBamFeaturesParams struct containing parameters for feature computation
   */
  void SetupWorkersAndAlignmentReaders(const ComputeBamFeaturesCliParams& cli_params,
                                       const ComputeBamFeaturesParams& params);

  /**
   * @brief Extracts BAM features in parallel for a list of genomic regions. Directly calls
   * ExtractBamFeatures to extract BAM features using the number of threads specified via the cli_params running in
   * parallel for the regions in the regions_list.
   * @param regions_list A list of genomic regions to process.
   * @param cli_params User-specified options for feature computation
   * @param progress Progress meter for logging and tracking
   */
  void ExecuteBamFeatureExtraction(const vec<Region>& regions_list,
                                   const ComputeBamFeaturesCliParams& cli_params,
                                   Progress& progress);

 private:
  /**
   * @brief Extracts BAM features for a given region using a specific thread and updates progress.
   *        This function is called in parallel for each region and is responsible for invoking
   *        the worker to compute features and updating the covered bases.
   * @param thread_id The id of the thread processing this region
   * @param region The genomic region to process
   * @param cli_params User-specified options for feature computation
   * @param progress Progress meter for logging and tracking
   */
  void ExtractBamFeatures(int thread_id,
                          const Region& region,
                          const ComputeBamFeaturesCliParams& cli_params,
                          Progress& progress);
  /**
   * @brief Breaks up a BAM into a list of regions. If bed regions are provided the bam is broken up based on the
   * regions specified, otherwise each chromosome is split based on the max_blocksize
   * @param contig_map A map of RefRegion structs containing the start and end aligned positions for each contig
   * @param bed_regions A list of bed regions corresponding to regions of interest
   * @param regions_queue A vector of regions that are of most max_region_size_per_thread in size on the bam
   * @param max_region_size_per_thread A maximum size each region can be
   */
  static void PopulateRegionsQueue(const StrMap<RefRegion>& contig_map,
                                   const StrMap<vec<Interval>>& bed_regions,
                                   vec<Region>& regions_queue,
                                   u64 max_region_size_per_thread);

  /**
   * @brief Breaks up a BAM into a list of regions corresponding to the given list of bed regions such that each regions
   * list of bam blocks is no larger than the max region size
   * @param contig_map A map of RefRegion structs containing the start and end aligned positions for each contig
   * @param bed_regions A list of bed regions corresponding to regions of interest
   * @param regions_queue A vector of regions that are of most max_region_size_per_thread in size on the bam. Each
   * region corresponds to at most one bed region. If a region is larger than max_region_size_per_thread, it is broken
   * up into multiple regions of size at most max_region_size_per_thread
   * @param max_region_size_per_thread A maximum size each region can be
   */
  static void GetGenomicRegionsForBedFile(const StrMap<RefRegion>& contig_map,
                                          const StrMap<vec<Interval>>& bed_regions,
                                          vec<Region>& regions_queue,
                                          u64 max_region_size_per_thread);

  /**
   * @brief Breaks up a BAM into a list of regions corresponding to groups of bam blocks that are no larger than the max
   * region size
   * @param contig_map A map of RefRegion structs containing the start and end aligned positions for each contig
   * @param regions_queue A vector of regions that are of most max_region_size_per_thread in size on the bam
   * @param max_region_size_per_thread A maximum size each region can be
   */
  static void GetGenomicRegionsForChromosome(const StrMap<RefRegion>& contig_map,
                                             vec<Region>& regions_queue,
                                             u64 max_region_size_per_thread);

  AlignmentReaderCache _alignment_reader_cache;
  vec<ComputeBamRegionFeatures> _workers;
  vec<vec<AlignmentReader>> _alignment_readers;
  StrMap<std::string> _ref_seqs{};
  std::atomic_uint32_t _bases_covered;
  std::unique_ptr<LockedTsvWriter> _writer;
};

}  // namespace xoos::svc
