#pragma once

#include <xoos/io/bed-region.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "compute-bam-features/alignment-reader.h"
#include "core/variant-info.h"
#include "ref-region.h"
#include "region.h"
#include "util/locked-tsv-writer.h"

namespace xoos::svc {

// invariant parameters for parallel feature extraction
struct ComputeBamFeaturesParams {
  vec<UnifiedFeatureCols> feature_cols;
  u8 min_bq{};
  u8 min_mapq{};
  float min_allowed_distance_from_end{};
  u32 min_family_size{};
  std::optional<u32> max_read_variant_count{};
  std::optional<float> max_read_variant_count_normalized{};
  u64 min_homopolymer_length{};
  bool duplex{};
  bool filter_homopolymer{};
  std::optional<std::string> tumor_read_group{};
  bool decode_yc{};
  yc_decode::BaseType min_base_type{};
};

using BamFeatures = std::tuple<UnifiedVariantFeatures, UnifiedReferenceFeatures>;

/**
 * @brief Worker class that computes variant and reference features from BAM alignments within a specific genomic
 * region. Processes individual reads and extracts features used for variant calling and filtering.
 */
class ComputeBamRegionFeatures {
 public:
  /**
   * @brief Constructs a ComputeBamRegionFeatures worker with the specified alignment readers, parameters, and output
   * writer.
   * @param alignment_readers Vector of alignment readers for accessing BAM files
   * @param params Parameters controlling feature computation behavior
   * @param writer Optional output writer for serializing computed features to file
   */
  ComputeBamRegionFeatures(const vec<AlignmentReader>& alignment_readers,
                           ComputeBamFeaturesParams params,
                           LockedTsvWriter* writer);

  /**
   * @brief Computes variant and reference feature sets for all reads aligning within the specified region.
   * @param region Positional data for the BAM block region
   * @param ref_sequence The reference chromosome sequence
   * @param chrom_to_vcf_features Map of chromosome to VCF features
   * @param skip_variants Set of known variants to skip
   * @return Variant features and reference allele features
   */
  BamFeatures ComputeBamFeatures(const Region& region,
                                 const std::string& ref_sequence,
                                 const std::optional<ChromToVcfFeaturesMap>& chrom_to_vcf_features,
                                 const std::optional<StrUnorderedSet>& skip_variants);

 public:
  // The following methods are made public to enable unit testing, but they are not intended
  // to be used directly.

 private:
  /**
   * @brief Map of read names to unique numeric identifiers for tracking reads across the region
   */
  using ReadNameToId = std::unordered_map<std::string, ReadId>;

  /**
   * @brief Processes and finalizes variant features after all reads have been analyzed.
   *        Computes summary statistics, handles multiple variants at the same position,
   *        and serializes features for output.
   */
  void ProcessVariantFeatures();

  /**
   * @brief Computes features for all alignments in a single BAM file within the specified region.
   * @param region The genomic region to process
   * @param ref_sequence The reference chromosome sequence
   * @param skip_variants Set of known variants to skip during processing
   * @param alignment_reader The BAM file reader to process
   */
  void ComputeFeaturesForBam(const Region& region,
                             const std::string& ref_sequence,
                             const std::optional<StrUnorderedSet>& skip_variants,
                             const AlignmentReader& alignment_reader);

  /**
   * @brief Computes features for a single read alignment within the specified region.
   *        Handles duplex vs non-duplex reads and YC tag decoding.
   * @param region The genomic region being processed
   * @param ref_sequence The reference chromosome sequence
   * @param skip_variants Set of known variants to skip during processing
   * @param b The BAM record to process
   */
  void ComputeFeaturesForRead(const Region& region,
                              const std::string& ref_sequence,
                              const std::optional<StrUnorderedSet>& skip_variants,
                              const bam1_t* b);

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

  /**
   * @brief Helper function to compute final feature values after all reads have been processed.
   *        Calculates allele frequencies, alignment bias, and tumor-normal specific features.
   * @param var_feat Variant feature structure to update
   * @param ref_feat Reference feature structure to update
   * @param total_mapq_sum Total mapping quality sum across all reads at this position
   * @param total_baseq_sum Total base quality sum across all reads at this position
   * @param tally_tumor_normal_features Whether to compute tumor-normal specific features
   */
  static void UpdateFeaturesHelper(UnifiedVariantFeature& var_feat,
                                   UnifiedReferenceFeature& ref_feat,
                                   u32 total_mapq_sum,
                                   double total_baseq_sum,
                                   bool tally_tumor_normal_features);

  /**
   * @brief Helper function to serialize feature values to strings for output.
   *        Converts variant and reference features to string representations based on selected columns.
   * @param vid Variant identifier containing position and allele information
   * @param var_feat Variant feature values to serialize
   * @param ref_feat Reference feature values to serialize
   * @param feature_cols Selected feature columns to include in output
   * @return Vector of feature strings in the order specified by feature_cols
   */
  static vec<std::string> SerializeFeatureHelper(const VariantId& vid,
                                                 const UnifiedVariantFeature& var_feat,
                                                 const UnifiedReferenceFeature& ref_feat,
                                                 const vec<UnifiedFeatureCols>& feature_cols);

  /**
   * @brief Helper function to compute the feature set for a single BAM record.
   * @param record The BAM record to compute features for
   * @param ref_seq The reference chromosome sequence
   * @param region The region of interest
   * @param read_names A map of read names to read ids
   * @param current_read_id The current read id to assign to the read name if it is not present in the map
   * @param skip_variants The set of variants to skip
   */
  void ComputeBamFeaturesHelper(const bam1_t* record,
                                const std::string& ref_seq,
                                const Region& region,
                                ReadNameToId& read_names,
                                ReadId& current_read_id,
                                const std::optional<StrUnorderedSet>& skip_variants);

  /**
   * @brief Get the read id for the given read name. If the read name is not present in the read_names map, it is added
   * to the map with the current_read_id as the value, current_read_id is then incremented.
   * @param [in/out] read_names A map of read names to read ids
   * @param [in] read_name The read name to get the id for
   * @param [in/out] current_read_id The current read id to assign to the read name if it is not present in the map
   * @return The read id for the given read name
   */
  static ReadId GetReadId(ReadNameToId& read_names, const std::string& read_name, ReadId& current_read_id);

  /**
   * @brief Resets all feature storage and tracking variables to prepare for processing a new genomic region.
   *        Clears read name mappings, variant/reference features, and VCF position caches.
   */
  void ResetForRegion();

 private:
  const vec<AlignmentReader>& _alignment_readers;
  ComputeBamFeaturesParams _params{};
  LockedTsvWriter* _writer{nullptr};

  // The following variables are used to track read alignments and their features,
  // they are reset for each genomic region.
  ReadNameToId _read_names{};
  ReadId _current_read_id{0};

  UnifiedVariantFeatures _var_features;
  // Store reference allele features at variant positions
  UnifiedReferenceFeatures _ref_features;
  PositionToVcfFeaturesMap _vcf_features;
  std::set<u64> _vcf_positions;
};

}  // namespace xoos::svc
