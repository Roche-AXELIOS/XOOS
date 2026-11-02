#include "intervals/augment-intervals.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <utility>

#include <fmt/core.h>

#include <xoos/io/fasta-reader.h>
#include <xoos/log/logging.h>
#include <xoos/sex_predict/sex.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-options.h"
#include "intervals/augmented-intervals.h"
#include "intervals/check-chr-y-par.h"
#include "intervals/interval-tools.h"
#include "intervals/intervals.h"
#include "io/bigwig.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/fai.h"
#include "segmentation/genomic-segments.h"
#include "segmentation/read-segments.h"
#include "segmentation/segment-type.h"

namespace xoos::cnc {

void AugmentIntervalsMain(const CopyNumberCallerOptions& options) {
  AugmentedIntervals augmented_intervals;
  const auto augmented_intervals_out = options.output_dir / kDefaultAugmentedIntervalsOutput;
  file::CheckOutputFilePathsForWriting({augmented_intervals_out}, options.overwrite);
  if (options.augment_intervals_options.generate_whole_genome_intervals) {
    augmented_intervals = GenerateAndAugmentWholeGenomeIntervals(options);
  } else {
    if (!options.augmented_intervals_fname.has_value()) {
      Logging::Error("Intervals file must be provided if not generating whole-genome intervals");
      throw std::runtime_error("Intervals file must be provided if not generating whole-genome intervals");
    }
    augmented_intervals = GenerateAndAugmentOnAndOffTargetIntervals(options);
  }
  const vec<GenomicSegment> seed_segments = ReadSegments(options.seed_segments_fname, SegmentType::kSeed);
  augmented_intervals = RemovePARIfChrYPARUnmasked(options.normal_bam_fname, augmented_intervals, seed_segments);
  std::ofstream ofs(augmented_intervals_out);
  augmented_intervals.Write(ofs, options.command_line_info);
}

AugmentedIntervals GenerateAndAugmentWholeGenomeIntervals(const CopyNumberCallerOptions& options) {
  SeqLenMap seq_lengths(ParseFai(options.reference_genome_fai_fname));

  BigWig mapp_bigwig =
      ValidateAndOpenMappabilityBigWig(options.mappability_bigwig_fname, options.reference_genome_fai_fname);
  ContigToIntervals intervals;
  f64 default_min_mappability = NAN;
  Logging::Info("Generating {}-bp intervals across whole genome",
                options.augment_intervals_options.whole_genome_interval_size);
  intervals = IntervalsFromSeqLenMap(seq_lengths, options.augment_intervals_options.include_alt_contigs);
  ThrowIfEmpty(intervals);
  Logging::Info("{} intervals", GetNContigToIntervals(intervals));
  Logging::Info("Removing blocklisted regions from intervals");
  if (options.blocklist_bed_fname.has_value()) {
    const ContigToIntervals blocklist_intervals = BedToMap(options.blocklist_bed_fname.value());
    intervals = SubtractRegionsFromIntervals(intervals, blocklist_intervals);
    Logging::Info("{} intervals after removing blocklisted regions", GetNContigToIntervals(intervals));
  }
  default_min_mappability = options.augment_intervals_options.whole_genome_min_mappability;
  Logging::Info("{} intervals", GetNContigToIntervals(intervals));
  Logging::Info("Cutting mappability=0 regions from intervals");
  intervals = Remove0MappableRegions(intervals, mapp_bigwig);
  ThrowIfEmpty(intervals, "No usable intervals. Check if mappability values overlap with reference genome");
  Logging::Info("{} intervals after removing 0-mappability regions", GetNContigToIntervals(intervals));
  Logging::Info("Splitting intervals into sub-intervals of length {}",
                options.augment_intervals_options.whole_genome_interval_size);
  intervals = SplitIntervalsByWindowSizedLastWindowSizeOfRemainder(
      intervals, options.augment_intervals_options.whole_genome_interval_size);
  ThrowIfEmpty(intervals);
  Logging::Info("{} intervals", GetNContigToIntervals(intervals));
  Logging::Info("Removing low mappability autosomal on-target intervals");
  intervals = RemoveLowMappableRegions(
      intervals, mapp_bigwig, default_min_mappability, sex_predict::SexChromHandling::kNoSexChrom);
  ThrowIfEmpty(intervals);
  Logging::Info("Removing low mappability allosomal on-target intervals");
  intervals = RemoveLowMappableRegions(intervals,
                                       mapp_bigwig,
                                       options.augment_intervals_options.sex_chromosome_min_mappability,
                                       sex_predict::SexChromHandling::kOnlySexChrom);
  ThrowIfEmpty(intervals);
  Logging::Info("{} intervals", GetNContigToIntervals(intervals));
  // now put intervals into an AugmentedIntervals class
  AugmentedIntervals augmented_intervals;
  augmented_intervals.SetRegions(ToRegionStringVec(intervals));
  augmented_intervals.SetOnTargetStatus(vec<bool>(GetNContigToIntervals(intervals), true));
  intervals = RemoveRegionsByLength(intervals, options.augment_intervals_options.whole_genome_interval_size);
  // Update augmented_intervals with modified intervals
  augmented_intervals.SetRegions(ToRegionStringVec(intervals));
  // Update on-target status
  augmented_intervals.SetOnTargetStatus(vec<bool>(GetNContigToIntervals(intervals), true));
  // annotate with original contig lengths
  augmented_intervals.SetSeqLengths(seq_lengths);
  Logging::Info("Sorting regions");
  // sort
  augmented_intervals.SortByRegion(options.reference_genome_fai_fname);
  io::FastaReader fa(options.reference_genome_fname);

  Logging::Info("Annotating intervals with GC bias");
  augmented_intervals.AnnotateGC(fa);
  Logging::Info("Removing intervals with extreme GC content");
  augmented_intervals.FilterExtremeGC(options.augment_intervals_options.min_gc_content,
                                      options.augment_intervals_options.max_gc_content);

  Logging::Info("Annotating intervals with mappability");
  augmented_intervals.AnnotateMappability(mapp_bigwig);
  augmented_intervals.SetReferenceFile(options.reference_genome_fname);
  return augmented_intervals;
}

AugmentedIntervals GenerateAndAugmentOnAndOffTargetIntervals(const CopyNumberCallerOptions& options) {
  SeqLenMap seq_lengths(ParseFai(options.reference_genome_fai_fname));
  BigWig mapp_bigwig =
      ValidateAndOpenMappabilityBigWig(options.mappability_bigwig_fname, options.reference_genome_fai_fname);
  ContigToIntervals on_targets;
  // if we're generating from whole genome, we also have to disable off-targets
  f64 default_min_mappability = NAN;
  Logging::Info("Generating on-targets from bed file");
  on_targets = BedToMap(options.augmented_intervals_fname.value());
  ThrowIfEmpty(on_targets);
  Logging::Info("{} intervals", GetNContigToIntervals(on_targets));
  ContigToIntervals blocklist_intervals;
  if (options.blocklist_bed_fname.has_value()) {
    Logging::Info("Removing blocklisted regions from on-target intervals");
    blocklist_intervals = BedToMap(options.blocklist_bed_fname.value());
    // NOLINTNEXTLINE(readability-suspicious-call-argument) - subtracting blocklist from targets; order is intentional
    on_targets = SubtractRegionsFromIntervals(on_targets, blocklist_intervals);
    Logging::Info("{} on-target intervals after removing blocklisted regions", GetNContigToIntervals(on_targets));
  }
  default_min_mappability = options.augment_intervals_options.on_target_min_mappability;
  Logging::Info("{} on-target intervals", GetNContigToIntervals(on_targets));
  Logging::Info("Cutting mappability=0 regions from on-targets");
  on_targets = Remove0MappableRegions(on_targets, mapp_bigwig);
  ThrowIfEmpty(on_targets);
  Logging::Info("{} on-target intervals after removing 0-mappability regions", GetNContigToIntervals(on_targets));
  Logging::Info("Splitting on-target intervals into sub-intervals of average length {}",
                options.augment_intervals_options.on_target_interval_size);
  on_targets = SplitIntervalsByWindowSizeAdjustForEqualWindows(
      on_targets, options.augment_intervals_options.on_target_interval_size);
  ThrowIfEmpty(on_targets);
  Logging::Info("{} on-target intervals", GetNContigToIntervals(on_targets));
  Logging::Info("Removing low mappability autosomal on-target intervals");
  on_targets = RemoveLowMappableRegions(
      on_targets, mapp_bigwig, default_min_mappability, sex_predict::SexChromHandling::kNoSexChrom);
  ThrowIfEmpty(on_targets);
  Logging::Info("Removing low mappability allosomal on-target intervals");
  on_targets = RemoveLowMappableRegions(on_targets,
                                        mapp_bigwig,
                                        options.augment_intervals_options.sex_chromosome_min_mappability,
                                        sex_predict::SexChromHandling::kOnlySexChrom);
  ThrowIfEmpty(on_targets);
  Logging::Info("{} on-target intervals", GetNContigToIntervals(on_targets));
  AugmentedIntervals augmented_intervals;
  if (!options.augment_intervals_options.no_off_targets) {
    Logging::Info("Generating off-target intervals");
    const ContigToIntervals off_targets =
        GenerateOffTargets(on_targets,
                           seq_lengths,
                           mapp_bigwig,
                           blocklist_intervals,
                           options.augment_intervals_options.off_target_min_width,
                           options.augment_intervals_options.off_target_interval_size,
                           options.augment_intervals_options.off_target_trim_length,
                           options.augment_intervals_options.off_target_min_mappability,
                           options.augment_intervals_options.sex_chromosome_min_mappability);
    // concat on- and off- target intervals. they will be sorted within the augmented_intervals object
    vec<std::string> on_target_regions = ToRegionStringVec(on_targets);
    vec<std::string> off_target_regions = ToRegionStringVec(off_targets);
    const size_t n_on_target = on_target_regions.size();
    const size_t n_off_target = off_target_regions.size();
    vec<std::string> all_regions;
    all_regions.reserve(n_on_target + n_off_target);
    std::ranges::move(on_target_regions, std::back_inserter(all_regions));
    std::ranges::move(off_target_regions, std::back_inserter(all_regions));
    // create the on-target std::vector
    vec<bool> on_target_status;
    on_target_status.reserve(all_regions.size());
    for (size_t i = 0; i < n_on_target; ++i) {
      on_target_status.push_back(true);
    }
    for (size_t i = 0; i < n_off_target; ++i) {
      on_target_status.push_back(false);
    }
    augmented_intervals.SetRegions(std::move(all_regions));
    augmented_intervals.SetOnTargetStatus(std::move(on_target_status));
  } else {
    Logging::Info("Off target generation was disabled");
    augmented_intervals.SetRegions(ToRegionStringVec(on_targets));
    augmented_intervals.SetOnTargetStatus(vec<bool>(GetNContigToIntervals(on_targets), true));
  }
  augmented_intervals.SetSeqLengths(seq_lengths);
  Logging::Info("Sorting regions");
  augmented_intervals.SortByRegion(options.reference_genome_fai_fname);
  io::FastaReader fa(options.reference_genome_fname);
  Logging::Info("Annotating intervals with GC bias");
  augmented_intervals.AnnotateGC(fa);
  Logging::Info("Annotating intervals with mappability");
  augmented_intervals.AnnotateMappability(mapp_bigwig);
  augmented_intervals.SetReferenceFile(options.reference_genome_fname);
  return augmented_intervals;
}
}  // namespace xoos::cnc
