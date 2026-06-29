#include "parallel-compute-utils.h"

#include <algorithm>

#include "compute-bam-features/region.h"
#include "compute-vcf-features/compute-vcf-features.h"
#include "vcf-to-bed/vcf-to-bed.h"

namespace xoos::svc {

WorkerContext::WorkerContext(const std::string& vcf_file,
                             const std::optional<std::string>& popaf_file,
                             const std::vector<fs::path>& bam_inputs,
                             AlignmentReaderCache& alignment_reader_cache)
    // Set up pointers to BAM readers, SAM headers, and BAM indexes for all input BAM files
    : filter_variants_reader{vcf_file}, vcf_feature_extraction_reader{vcf_file} {
  alignment_readers = alignment_reader_cache.Open(bam_inputs, "r");
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
vec<TargetRegion> PartitionVcfRegions(const fs::path& vcf_file,
                                      const size_t threads,
                                      std::optional<ChromIntervalsMap> bed_regions) {
  vec<TargetRegion> partitioned_regions;
  const u32 left_pad = 1;
  const u32 right_pad = 1;
  ChromIntervalsMap target_regions;
  if (bed_regions.has_value()) {
    target_regions = bed_regions.value();
  }

  // Only collapse intervals if they are immediately adjacent to each other. So, each chromosome would have
  // approximately one interval for each variant. The only exception is large deletions, which have large intervals.
  const u32 collapse_dist = 0;

  // Use a single thread to extract variant positions. If the VCF file has many chromosomes, then the parallelized
  // version (1 parallel task per chromomosome) may cause this step to become very slow.
  auto chrom_to_intervals =
      ExtractVariantIntervalsSingleThreaded(vcf_file, target_regions, left_pad, right_pad, collapse_dist);
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

/**
 * @brief Computes BAM and VCF features for a given region.
 * @param global_ctx Global context containing configuration and reference sequences.
 * @param worker_ctx Worker context containing BAM readers, SAM headers, and BAM indexes.
 * @param region The target region for which to compute features.
 * @param bed_regions Optional map of chromosome to vector of intervals representing BED regions.
 * @param interest_regions Optional map of chromosome to vector of intervals representing regions of interest.
 * @return A tuple containing the computed VCF features and BAM features.
 */
std::tuple<VarIdToVcfFeatures, BamRegionFeatureCollection> ComputeBamAndVcfFeaturesForRegion(
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
  VarIdToVcfFeatures vcf_features;
  BamRegionFeatureCollection bam_features;
  ChromPosToRefBamFeatures ref_features;
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
    std::optional<VarIdToVcfFeatures> passed_vcf_features = std::nullopt;
    if (use_vcf_features) {
      vcf_features = ExtractFeaturesForRegion(worker_ctx.vcf_feature_extraction_reader,
                                              global_ctx.genome,
                                              worker_ctx.popaf_reader,
                                              target_regions,
                                              interest_regions,
                                              region.chrom,
                                              global_ctx.is_germline_tagging);
      passed_vcf_features = vcf_features;
    }

    bam_features =
        ComputeBamRegionFeatures(worker_ctx.alignment_readers, global_ctx.bam_feat_params, nullptr)
            .ComputeBamFeatures(bam_region, global_ctx.ref_seqs.at(region.chrom), passed_vcf_features, skip_variants);
  }
  return std::make_tuple(vcf_features, bam_features);
}

}  // namespace xoos::svc
