#include "check-chr-y-par.h"

#include <algorithm>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "util/eigen-math.h"
#include "util/parse-util.h"

namespace xoos::cnc {
/**
 * @brief Checks whether the chrY PAR region was masked during alignment by looking for any alignments in the chrY PAR
 * region defined in the seed segments file. If any alignments are found in the chrY PAR region, an error is thrown
 * since this could significantly affect copy number predictions in allosomes. If the seed segments do not specify any
 * segments in the chrY PAR region, a warning is logged since we cannot check whether the chrY PAR region was masked
 * during alignment.
 */
bool CheckChrYPARUnmasked(const fs::path& bam_file, const vec<GenomicSegment>& seed_segments) {
  // find the chrY seed segment that is in PAR region
  const auto result = std::ranges::find_if(seed_segments, [](const GenomicSegment& seg) {
    return ((seg.contig == "chrY" || seg.contig == "Y") && seg.in_pseudo_autosomal_region.has_value() &&
            seg.in_pseudo_autosomal_region.value());
  });
  // no need to throw error if there aren't any seed segments in chrY PAR region, as that likely means the user is not
  // using any seed segments in the chrY PAR region for calling (which is reasonable since those regions can be noisy)
  if (result == seed_segments.end()) {
    Logging::Warn(
        "Could not find any seed segments in chrY PAR region. Cannot check whether chrY PAR region was masked during "
        "alignment. This could affect copy number predictions in allosomes.");
  } else {
    // open the BAM file and check for reads that are in the chrY PAR region described by the seed segment file.
    const io::SamFilePtr sam_file = io::SamOpen(bam_file.string(), "r");
    const io::SamHdrPtr sam_hdr = io::SamHdrRead(sam_file.get());
    const fs::path bai_file = bam_file.string() + ".bai";
    const io::HtsIdxPtr sam_idx = io::SamIndexLoad(sam_file.get(), bai_file);
    const std::string region_string =
        result->contig + ":" + std::to_string(result->start + 1) + "-" + std::to_string(result->end);
    const io::HtsItrPtr sam_itr_ptr = io::SamItrQuerySNoThrow(sam_idx.get(), sam_hdr.get(), region_string);
    // if null is returned by sam_itr_query, then it means that chrY or the PAR region was not in the reference FASTA
    // used by the aligner. Not issuing a warning because it is not deletrious behavior.
    if (sam_itr_ptr.get() == nullptr) {
      Logging::Info("BAM file has no alignments in the chrY PAR regions");
    } else {
      size_t n_alignments = 0;
      const io::Bam1Ptr bam1_ptr(bam_init1());
      while (sam_itr_next(sam_file.get(), sam_itr_ptr.get(), bam1_ptr.get()) > 0) {
        // only increment n_alignments if the read starts AFTER the start of the PAR region and ends BEFORE the end of
        // the PAR region, since there can be alignments that partially overlap the PAR region but are not fully
        // contained within it
        if (bam1_ptr->core.pos >= ToSigned(result->start) && bam_endpos(bam1_ptr.get()) <= ToSigned(result->end)) {
          ++n_alignments;
        }
      }
      if (n_alignments > 0) {
        Logging::Warn(
            "Found {} alignments in chrY PAR region. This suggests that the chrY PAR region may NOT have been masked "
            "during alignment. This could affect copy number predictions in allosomes.",
            n_alignments);
        return true;
      }
    }
  }
  return false;
}

namespace {
AugmentedIntervals RemovePARFromIntervals(const AugmentedIntervals& augmented_intervals,
                                          const vec<GenomicSegment>& seed_segments) {
  // Get the PAR regions from the seed segments
  vec<std::tuple<std::string, size_t, size_t>> par_regions;
  for (const auto& seg : seed_segments) {
    if (seg.in_pseudo_autosomal_region.has_value() && seg.in_pseudo_autosomal_region.value()) {
      par_regions.emplace_back(seg.contig, seg.start, seg.end);
    }
  }
  if (par_regions.empty()) {
    Logging::Warn("No seed segments found in chrY PAR region. Cannot remove PAR regions from augmented intervals.");
    return augmented_intervals;
  }
  // Remove intervals that overlap with the PAR regions
  vec<std::string> filtered_regions;
  vec<f64> filtered_gc_bias;
  vec<f64> filtered_mappability;
  vec<bool> filtered_on_target;
  const auto& regions = augmented_intervals.GetRegions();
  const auto& gc_biases = augmented_intervals.GetGCBias();
  const auto& mappabilities = augmented_intervals.GetMappability();
  const auto& on_targets = augmented_intervals.GetOnTargetStatus();
  for (size_t i = 0; i < augmented_intervals.GetRegions().size(); ++i) {
    // parse the region string to get the contig, start, and end
    auto [contig, start, end] = ParseRegionString(regions[i]);
    const bool in_par_region = std::ranges::any_of(
        par_regions, [&contig, start, end](const std::tuple<std::string, size_t, size_t>& par_region) {
          const auto& [par_contig, par_start, par_end] = par_region;
          return (contig == par_contig && (start <= par_end && end >= par_start));
        });
    if (!in_par_region) {
      filtered_regions.push_back(regions[i]);
      filtered_gc_bias.push_back(gc_biases[static_cast<Eigen::Index>(i)]);
      filtered_mappability.push_back(mappabilities[static_cast<Eigen::Index>(i)]);
      filtered_on_target.push_back(on_targets[i]);
    }
  }
  AugmentedIntervals filtered_augmented_intervals;
  filtered_augmented_intervals.SetRegions(std::move(filtered_regions));
  filtered_augmented_intervals.SetGCBias(eigen::FromStdVec(filtered_gc_bias));
  filtered_augmented_intervals.SetMappability(eigen::FromStdVec(filtered_mappability));
  filtered_augmented_intervals.SetOnTargetStatus(std::move(filtered_on_target));
  auto seq_lengths = augmented_intervals.GetSeqLengths();
  filtered_augmented_intervals.SetSeqLengths(seq_lengths);
  filtered_augmented_intervals.SetReferenceFile(augmented_intervals.GetReferenceFile());
  Logging::Info(
      "Removed intervals that overlap with PAR regions. Original number of intervals: {}, number of intervals after "
      "filtering: {}.",
      augmented_intervals.GetRegions().size(),
      filtered_augmented_intervals.GetRegions().size());
  return filtered_augmented_intervals;
}
}  // namespace

AugmentedIntervals RemovePARIfChrYPARUnmasked(const fs::path& bam_file,
                                              const AugmentedIntervals& augmented_intervals,
                                              const vec<GenomicSegment>& seed_segments) {
  if (CheckChrYPARUnmasked(bam_file, seed_segments)) {
    return RemovePARFromIntervals(augmented_intervals, seed_segments);
  }
  Logging::Info("chrY PAR region appears to be masked during alignment. PAR regions will be included in analysis");
  return augmented_intervals;
}
}  // namespace xoos::cnc
