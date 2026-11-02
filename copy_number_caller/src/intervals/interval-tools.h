#pragma once
#include <optional>
#include <string>

#include <xoos/io/bed-region.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

#include "intervals/intervals.h"
#include "io/bigwig.h"
#include "io/fai.h"

namespace xoos::cnc {

ContigToIntervals MergeAdjacentIntervals(const ContigToIntervals& intervals_map);
ContigToIntervals GenerateCandidateOffTargets(const ContigToIntervals& intervals_map,
                                              const StrUnorderedMap<size_t>& seq_lens);
ContigToIntervals RemoveLowMappableRegions(
    const ContigToIntervals& map,
    BigWig& mapp_bigwig,
    f64 min_mappability,
    sex_predict::SexChromHandling sex_chrom_behavior = sex_predict::SexChromHandling::kNoSexChrom);
ContigToIntervals RemoveRegionsByLength(const ContigToIntervals& map, size_t length);
ContigToIntervals SplitIntervalsByWindowSizeAdjustForEqualWindows(const ContigToIntervals& map, size_t interval_size);
ContigToIntervals BedToMap(const fs::path& bed_file);
ContigToIntervals GenerateOffTargets(const ContigToIntervals& on_targets,
                                     const SeqLenMap& seq_lengths,
                                     BigWig& mapp_bigwig,
                                     const ContigToIntervals& blocklist_intervals,
                                     size_t off_target_min_width,
                                     size_t off_target_interval_size,
                                     size_t off_target_trim_length,
                                     f64 off_target_min_mappability,
                                     f64 sex_target_min_mappability);
size_t GetNContigToIntervals(const ContigToIntervals& map);
vec<std::string> ToRegionStringVec(const ContigToIntervals& map);
ContigToIntervals Remove0MappableRegions(const ContigToIntervals& map, BigWig& mapp_bigwig);
ContigToIntervals SplitIntervalsByWindowSizedLastWindowSizeOfRemainder(const ContigToIntervals& map,
                                                                       size_t interval_size);
ContigToIntervals IntervalsFromSeqLenMap(const SeqLenMap& seq_lengths, bool include_alt_contigs);
vec<Interval> SubtractIntervalBFromIntervalA(const Interval& a, const Interval& b);
ContigToIntervals SubtractRegionsFromIntervals(const ContigToIntervals& intervals,
                                               const ContigToIntervals& intervals_to_remove);
void ThrowIfEmpty(const ContigToIntervals& map, const std::optional<std::string>& diagnostic_message = std::nullopt);
}  // namespace xoos::cnc
