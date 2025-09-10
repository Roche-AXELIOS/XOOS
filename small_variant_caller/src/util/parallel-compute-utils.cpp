#include "parallel-compute-utils.h"

#include <algorithm>

#include <htslib/bgzf.h>

#include <xoos/error/error.h>

#include "compute-bam-features/region.h"
#include "compute-vcf-features/compute-vcf-features.h"
#include "core/variant-feature-extraction.h"
#include "vcf-to-bed/vcf-to-bed.h"

namespace xoos::svc {

WorkerContext::WorkerContext(const std::string& vcf_file,
                             const std::optional<std::string>& popaf_file,
                             const std::vector<fs::path>& bam_inputs,
                             const SVCConfig& model_config,
                             AlignmentReaderCache& alignment_reader_cache)
    // Set up pointers to BAM readers, SAM headers, and BAM indexes for all input BAM files
    : filter_variants_reader{vcf_file}, vcf_feature_extraction_reader{vcf_file} {
  alignment_readers = alignment_reader_cache.Open(bam_inputs, "r");
  // Set up workflow-specific model calculators
  switch (model_config.workflow) {
    case Workflow::kGermlineMultiSample:
    case Workflow::kGermline: {
      calculators.emplace_back(model_config.snv_model_file, GetFeatureVecLength(model_config.snv_scoring_cols));
      calculators.emplace_back(model_config.indel_model_file, GetFeatureVecLength(model_config.indel_scoring_cols));
      break;
    }
#ifdef SOMATIC_ENABLE
    case Workflow::kSomatic: {
      calculators.emplace_back(model_config.model_file, GetFeatureVecLength(model_config.scoring_cols));
      break;
    }
    case Workflow::kSomaticTumorNormal: {
      calculators.emplace_back(model_config.model_file, GetFeatureVecLength(model_config.scoring_cols));
      calculators.emplace_back(model_config.germline_fail_model_file,
                               GetFeatureVecLength(model_config.germline_fail_scoring_cols));
      break;
    }
#endif  // SOMATIC_ENABLE
    default: {
      // Workflow not supported
      break;
    }
  }
  // Set up the population allele frequency VCF reader if provided
  if (popaf_file.has_value()) {
    popaf_reader = io::VcfReader(popaf_file.value());
  }
}

/**
 * @brief Partition a VCF file into regions based on the number of threads available.
 * Each region will contain approximately the same number of variants, allowing for parallel processing.
 * @param vcf_file Path to the VCF file to partition.
 * @param threads Number of threads to use for parallel processing.
 * @return A vector of TargetRegion objects representing the partitioned regions.
 */
vec<TargetRegion> PartitionVcfRegions(const fs::path& vcf_file, const size_t threads) {
  vec<TargetRegion> partitioned_regions;
  const u32 left_pad = 1;
  const u32 right_pad = 1;

  // Only collapse intervals if they are immediately adjacent to each other. So, each chromosome would have
  // approximately one interval for each variant. The only exception is large deletions, which have large intervals.
  const u32 collapse_dist = 0;

  // Use a single thread to extract variant positions. If the VCF file has many chromosomes, then the parallelized
  // version (1 parallel task per chromomosome) may cause this step to become very slow.
  auto chrom_to_intervals = ExtractVariantIntervalsSingleThreaded(vcf_file, {}, left_pad, right_pad, collapse_dist);
  // The goal here is to assign roughly same number of variants per parallel task in each chromosome. Since features
  // have already been computed, there is no need to create parallel tasks for small equal-sized regions across the
  // entire genome. Otherwise, many tasks would be created for regions that do not overlap any variants at all.
  for (auto& [chrom, intervals] : chrom_to_intervals) {
    if (!intervals.empty()) {
      auto num_intervals = intervals.size();
      if (num_intervals > threads) {
        // All partitions of this chromosome would contain near-equal number of intervals.
        auto rem = num_intervals % threads;
        if (rem != 0) {
          rem = 1;
        }
        auto num_intervals_per_thread = num_intervals / threads + rem;
        auto last_i = num_intervals - 1;
        for (size_t i = 0; i < num_intervals; i += num_intervals_per_thread) {
          auto max_i = std::min(last_i, i + num_intervals_per_thread - 1);
          partitioned_regions.emplace_back(chrom, intervals[i].start, intervals[max_i].end);
        }
      } else {
        // Too few intervals in this chromosome.
        // Create only one partition for all intervals in this chromosome.
        partitioned_regions.emplace_back(chrom, intervals[0].start, (intervals.end() - 1)->end);
      }
    }
  }
  return partitioned_regions;
}

#ifdef SOMATIC_ENABLE
/**
 * @brief Breaks up a BAM into a list of regions corresponding to groups of bam blocks that are no larger than the max
 * region size
 * @param contig_map A map of RefRegion structs containing the start and end aligned positions for each contig
 * @param partitioned_regions A vector of regions that are of most max_region_size_per_thread in size on the bam
 * @param max_region_size_per_thread A maximum size each region can be
 */
static void GetGenomicRegionsForChromosome(const StrMap<RefRegion>& contig_map,
                                           vec<TargetRegion>& partitioned_regions,
                                           const u64 max_region_size_per_thread) {
  // This function mirrors logic that is currently implemented in compute-bam-features
  // Can be removed once need to have somatic region splitting match compute-bam-features is no longer there
  // TODO : Remove and unify with germline region splitting once max_variants_in_read filter is fixed
  for (const auto& [chrom, region] : contig_map) {
    // leftmost aligned position of the contig
    const u64 region_start = region.start_position;
    // rightmost aligned position of the contig
    const u64 region_end = region.end_position;
    for (u64 start = region_start; start <= region_end + 1; start += max_region_size_per_thread) {
      const u64 end = std::min(start + max_region_size_per_thread, region_end + 1);
      partitioned_regions.emplace_back(chrom, start, end);
    }
  }
}

/**
 * @brief Breaks up a BAM into a list of regions corresponding to the given list of bed regions such that each regions
 * list of bam blocks is no larger than the max region size
 * @param contig_map A map of RefRegion structs containing the start and end aligned positions for each contig
 * @param bed_regions A list of bed regions corresponding to regions of interest
 * @param partitioned_regions A vector of regions that are of most max_region_size_per_thread in size on the bam. Each
 * region corresponds to at most one bed region. If a region is larger than max_region_size_per_thread, it is broken up
 * into multiple regions of size at most max_region_size_per_thread
 * @param max_region_size_per_thread A maximum size each region can be
 */
static void GetGenomicRegionsForBedFile(const StrMap<RefRegion>& contig_map,
                                        const StrMap<vec<Interval>>& bed_regions,
                                        vec<TargetRegion>& partitioned_regions,
                                        const u64 max_region_size_per_thread) {
  // This function mirrors logic that is currently implemented in compute-bam-features
  // Can be removed once need to have somatic region splitting match compute-bam-features is no longer there
  // TODO : Remove and unify with germline region splitting once max_variants_in_read filter is fixed
  for (const auto& [chrom, intervals] : bed_regions) {
    // skip if the contig has no alignments
    if (!contig_map.contains(chrom)) {
      continue;
    }
    for (const auto& interval : intervals) {
      const u64 region_start = interval.start;
      const u64 region_end = interval.end;
      for (u64 start = region_start; start <= region_end + 1; start += max_region_size_per_thread) {
        const u64 end = std::min(start + max_region_size_per_thread, region_end + 1);
        partitioned_regions.emplace_back(chrom, start, end);
      }
    }
  }
}

/**
 * @brief Partition regions for somatic workflow based on the BAM file and optional BED regions.
 * @param region_size The maximum size of each region to partition.
 * @param bam_file Pointer to the opened BAM file.
 * @param header Pointer to the SAM header of the BAM file.
 * @param idx Pointer to the BAM index.
 * @param bed_regions Optional map of chromosome to vector of intervals representing BED regions.
 * @return A vector of TargetRegion objects representing the partitioned regions.
 */
vec<TargetRegion> PartitionRegionsForSomatic(const u64 region_size,
                                             const HtsFilePtr& bam_file,
                                             const SamHdrPtr& header,
                                             const HtsIdxPtr& idx,
                                             const StrMap<vec<Interval>>& bed_regions) {
  // This function mirrors logic that is currently implemented in compute-bam-features
  // Can be updated once need to have somatic region splitting match compute-bam-features is no longer there
  // TODO: Update how forcecalls are done so we can limit the space where we need VCF and BAM features
  // TODO: unify with germline region splitting once max_variants_in_read filter is fixed
  StrMap<RefRegion> contig_map = GetRefRegions(bam_file, header, idx);
  vec<TargetRegion> partitioned_regions;
  if (!bed_regions.empty()) {
    GetGenomicRegionsForBedFile(contig_map, bed_regions, partitioned_regions, region_size);
  } else {
    GetGenomicRegionsForChromosome(contig_map, partitioned_regions, region_size);
  }
  return partitioned_regions;
}
#endif  // SOMATIC_ENABLE

/**
 * @brief Computes BAM and VCF features for a given region.
 * @param global_ctx Global context containing configuration and reference sequences.
 * @param worker_ctx Worker context containing BAM readers, SAM headers, and BAM indexes.
 * @param region The target region for which to compute features.
 * @param bed_regions Optional map of chromosome to vector of intervals representing BED regions.
 * @param interest_regions Optional map of chromosome to vector of intervals representing regions of interest.
 * @return A tuple containing the computed VCF features, BAM features, and reference features.
 */
std::tuple<ChromToVcfFeaturesMap, ChromToVariantInfoMap, RefInfoMap> ComputeBamAndVcfFeaturesForRegion(
    const GlobalContext& global_ctx,
    WorkerContext& worker_ctx,
    const TargetRegion& region,
    const ChromIntervalsMap& bed_regions,
    const ChromIntervalsMap& interest_regions) {
  // This function is designed to be used by an individual worker thread as part of a larger taskflow workflow where
  // each task is getting features for a given region.

  // Compute VCF features first. Once VCF features have been computed we can use the feature set when computing BAM
  // features. This should reduce peak memory usage and lower compute time.
  // If a set of BED Regions has been passed use the intersection of the target regions and the BAM region
  auto bam_region = Region{.chrom = region.chrom, .start = region.start, .end = region.end};
  ChromIntervalsMap target_regions;
  auto i_region = Interval{.start = region.start, .end = region.end};
  if (bed_regions.contains(region.chrom)) {
    for (const auto& interval : bed_regions.at(region.chrom)) {
      if (IntervalOverlap(interval, i_region)) {
        // Adjust interval if it falls partly outside the region
        auto interval_to_use = Interval{.start = interval.start, .end = interval.end};
        if (interval_to_use.start < i_region.start) {
          interval_to_use.start = i_region.start;
        }
        if (interval_to_use.end > i_region.end) {
          interval_to_use.end = i_region.end;
        }
        target_regions[region.chrom].emplace_back(interval_to_use);
      }
    }
  } else if (bed_regions.empty()) {
    target_regions[region.chrom].emplace_back(i_region);
  }
  ChromToVcfFeaturesMap vcf_features;
  std::optional<ChromToVcfFeaturesMap> passed_vcf_features = std::nullopt;
  ChromToVariantInfoMap bam_features;
  RefInfoMap ref_features;
  bool use_vcf_features = true;
  if (!global_ctx.model_config.HasVcfFeatureScoringCols() || !global_ctx.model_config.use_vcf_features) {
    use_vcf_features = false;
  }
  if (!target_regions.empty()) {
    std::optional<StrUnorderedSet> skip_variants = std::nullopt;
    if (global_ctx.skip_variants_vcf.has_value()) {
      skip_variants = std::make_optional(ExtractVariantKeySet(global_ctx.skip_variants_vcf.value(), bam_region));
    }

    // Only compute VCF features for the region if there is a target entry and we have to based on workflow and scoring
    // cols
    if (use_vcf_features) {
      vcf_features[region.chrom] = ExtractFeaturesForRegion(worker_ctx.vcf_feature_extraction_reader,
                                                            global_ctx.genome,
                                                            worker_ctx.popaf_reader,
                                                            target_regions,
                                                            interest_regions,
                                                            region.chrom);
      passed_vcf_features = vcf_features;
    }

    auto [vid_to_var_feat, pos_to_ref_feat] =
        ComputeBamRegionFeatures(worker_ctx.alignment_readers, global_ctx.bam_feat_params, nullptr)
            .ComputeBamFeatures(bam_region, global_ctx.ref_seqs.at(region.chrom), passed_vcf_features, skip_variants);
    for (const auto& [vid, var_feat] : vid_to_var_feat) {
      bam_features[vid.chrom][vid.pos][vid] = var_feat;
    }
    ref_features[region.chrom] = pos_to_ref_feat;
  }
  return std::make_tuple(vcf_features, bam_features, ref_features);
}

}  // namespace xoos::svc
