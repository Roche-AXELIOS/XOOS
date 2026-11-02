#include "io/bed-reader.h"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>

#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/io/bed-region.h>
#include <xoos/io/fast-bed-parser.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

#include "core/region-labels-set.h"

namespace xoos::alignment_metrics {

namespace {

/**
 * Parse the comma-separated label list in `label_field`, register each label in `registry`, and
 * return the OR of their single-bit masks. Empty labels between commas are rejected.
 */
u16 ParseLabels(const std::string_view label_field, RegionLabelsSet& registry) {
  u16 labels = 0;
  size_t pos = 0;
  while (pos <= label_field.size()) {
    const size_t comma = label_field.find(',', pos);
    const size_t end = comma == std::string_view::npos ? label_field.size() : comma;
    const std::string_view token = label_field.substr(pos, end - pos);
    labels = static_cast<u16>(labels | registry.AddLabel(token));
    if (comma == std::string_view::npos) {
      break;
    }
    pos = comma + 1;
  }
  return labels;
}

/**
 * Check if the given chromosome is not present in the BAM file header.
 * Returns true if the chromosome is unknown, false otherwise.
 * Throws if the BAM header cannot be parsed.
 */
bool IsUnknownContig(const std::string& chromosome, sam_hdr_t* sam_hdr_ptr) {
  const auto tid = sam_hdr_name2tid(sam_hdr_ptr, chromosome.c_str());
  constexpr s32 kTidNotFound = -1;
  constexpr s32 kTidHeaderParseError = -2;
  if (tid == kTidNotFound) {
    return true;
  }
  if (tid == kTidHeaderParseError) {
    throw error::Error("Could not parse header of BAM");
  }
  if (tid < 0) {
    throw error::Error("Could not determine tid of '{}' in BAM", chromosome);
  }
  return false;
}

}  // namespace

/**
 * @brief Reads genomic regions from a BED file and validates them against a BAM file header.
 *
 * Parses a BED file to extract genomic intervals (chromosome, start, end) and validates
 * that each chromosome exists in the BAM file header. Regions on unknown chromosomes are skipped with
 * a warning. BED-format validation (column count, numeric non-negative coordinates, start <= end) and
 * skipping of UCSC browser configuration lines are handled by io::ForEachBedRegionInFileFast. Validating
 * against the BAM header prevents processing regions that don't exist in the sequencing data, which
 * would waste computation and could cause errors. The returned regions are sorted by position within
 * each chromosome for efficient downstream processing. Empty BED files yield an empty map.
 *
 * @param bed_filename Path to the BED file containing regions of interest
 * @param bam_filename Path to the BAM file used for chromosome validation
 * @param registry Region labels set populated with labels encountered when `stratify_by_label` is set
 * @param stratify_by_label When true, parse the 4th column as a comma-separated list of region labels
 * @return Map of chromosome names to sorted vectors of genomic intervals
 * @throws error::Error if the BED file doesn't exist or contains invalid data
 */
RegionsByChromosome GetRegionsFromBed(const fs::path& bed_filename,
                                      const fs::path& bam_filename,
                                      RegionLabelsSet& registry,
                                      const bool stratify_by_label) {
  if (!exists(bed_filename)) {
    throw error::Error("BED file '{}' does not exist", bed_filename.string());
  }
  const auto hts_file_ptr = io::HtsOpen(bam_filename, "rb");
  const auto sam_hdr_ptr = io::SamHdrRead(hts_file_ptr.get());

  RegionsByChromosome regions;
  StrUnorderedSet unknown_contigs;
  u32 region_count = 0;

  auto add_bed_region_to_map = [&unknown_contigs, &sam_hdr_ptr, stratify_by_label, &registry, &regions, &region_count](
                                   const io::BedRegion& region) {
    if (unknown_contigs.contains(region.chromosome)) {
      return;
    }
    // Validate that the chromosome specified exists in the BAM header
    if (IsUnknownContig(region.chromosome, sam_hdr_ptr.get())) {
      Logging::Warn("Ignoring chromosome '{}' from BED, not found in BAM", region.chromosome);
      unknown_contigs.emplace(region.chromosome);
      return;
    }
    u16 labels = 0;
    if (stratify_by_label) {
      // The label column is captured by ForEachBedRegionInFileFast as the first extra column (column 4).
      if (region.data.empty()) {
        throw error::Error(
            "--stratify-by-label requires a 4th column with region labels, but region '{}:{}-{}' in the provided "
            "BED has none",
            region.chromosome,
            region.start,
            region.end);
      }
      labels = ParseLabels(region.data.front(), registry);
    }
    regions[region.chromosome].emplace_back(region.start, region.end, labels);
    ++region_count;
  };

  io::ForEachBedRegionInFileFast(bed_filename, add_bed_region_to_map);

  // Sort the regions for each chromosome
  for (auto& [chromosome, intervals] : regions) {
    if (!std::ranges::is_sorted(intervals)) {
      std::ranges::sort(intervals);
    }
  }
  Logging::Info("Read {} regions from BED file '{}'", region_count, bed_filename.string());
  return regions;
}

RegionsByChromosome GetRegionsFromBed(const fs::path& bed_filename, const fs::path& bam_filename) {
  RegionLabelsSet registry;
  return GetRegionsFromBed(bed_filename, bam_filename, registry, false);
}

/**
 * @brief Extracts all reference contigs from a BAM file header as full-length genomic regions.
 *
 * Creates regions spanning the entire length of each chromosome/contig defined in the BAM header.
 * This is used when no BED file is provided, allowing metrics to be calculated across all sequenced
 * regions. Each region starts at position 0 and extends to the full contig length.
 *
 * @param bam_filename Path to the BAM file to extract contig information from
 * @return Map of chromosome names to full-length genomic intervals
 */
RegionsByChromosome GetRegionsFromBam(const fs::path& bam_filename) {
  const auto hts_file_ptr = io::HtsOpen(bam_filename, "rb");
  const auto sam_hdr_ptr = io::SamHdrRead(hts_file_ptr.get());

  RegionsByChromosome regions;
  for (s32 i = 0; i < sam_hdr_ptr->n_targets; ++i) {
    const std::string chromosome(sam_hdr_tid2name(sam_hdr_ptr.get(), i));
    const u64 len = ToUnsigned(sam_hdr_tid2len(sam_hdr_ptr.get(), i));
    regions[chromosome].emplace_back(0, len);
  }

  Logging::Info("Read {} contigs from BAM file '{}'", regions.size(), bam_filename.string());
  return regions;
}

}  // namespace xoos::alignment_metrics
