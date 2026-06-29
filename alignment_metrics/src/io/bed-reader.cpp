#include "io/bed-reader.h"

#include <filesystem>
#include <string>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>

namespace xoos::alignment_metrics {

/**
 * @brief Reads genomic regions from a BED file and validates them against a BAM file header.
 *
 * Parses a BED file to extract genomic intervals (chromosome, start, end) and validates
 * that each chromosome exists in the BAM file header. Regions on unknown chromosomes are skipped with
 * a warning. The function performs validation to ensure BED format correctness (at least 3 columns,
 * numeric coordinates, non-negative values, start < end) and handles special UCSC browser configuration
 * lines. This validation is essential to prevent processing regions that don't exist in the sequencing
 * data, which would waste computation and could cause errors. The function returns regions sorted by
 * position within each chromosome for efficient downstream processing. Empty BED files are handled
 * gracefully by returning an empty map.
 *
 * @param bed_filename Path to the BED file containing regions of interest
 * @param bam_filename Path to the BAM file used for chromosome validation
 * @return Map of chromosome names to sorted vectors of genomic intervals
 * @throws error::Error if the BED file doesn't exist or contains invalid data
 */
RegionsByChromosome GetRegionsFromBed(const fs::path& bed_filename, const fs::path& bam_filename) {
  if (!exists(bed_filename)) {
    throw error::Error("BED file '{}' does not exist", bed_filename.string());
  }
  const auto hts_file_ptr = io::HtsOpen(bam_filename, "rb");
  const auto sam_hdr_ptr = io::SamHdrRead(hts_file_ptr.get());

  RegionsByChromosome regions;
  StrUnorderedSet unknown_chromosomes;

  // Set up the TSV read for reading the BED file
  const csv::CSVFormat format = csv::CSVFormat().delimiter('\t').no_header().trim({' '});
  if (fs::file_size(bed_filename) == 0) {
    return {};
  }
  csv::CSVReader reader(bed_filename.string(), format);
  u32 region_count = 0;
  for (const auto& row : reader) {
    auto chrom = row[0].get<std::string>();
    // skip two special cases which are used to configure UCSC genome browser
    if (chrom.starts_with("browser") || chrom.starts_with("track") || chrom.starts_with("#")) {
      continue;
    }
    if (row.size() < 3) {
      throw std::runtime_error("BED file must have at least 3 columns");
    }
    if (row[1].is_str() || row[2].is_str()) {
      throw std::runtime_error("BED file must have numeric start and end columns");
    }
    if (row[1].get<s64>() < 0 || row[2].get<s64>() < 0) {
      throw std::runtime_error("BED file must have non-negative start and end columns");
    }
    if (row[1].get<s64>() >= row[2].get<s64>()) {
      throw std::runtime_error("BED file must have start < end");
    }
    const auto start = row[1].get<s64>();
    const auto end = row[2].get<s64>();
    if (unknown_chromosomes.contains(chrom)) {
      continue;
    }
    // Validate that the chromosome specified exists in the BAM header
    const auto targeted_id = sam_hdr_name2tid(sam_hdr_ptr.get(), chrom.c_str());
    if (targeted_id == -1) {
      // Log a warning and skip this chromosome, add it to the unknown set to avoid logging many times
      Logging::Warn("Ignoring chromosome '{}' from '{}', not found in '{}'", chrom, bed_filename, bam_filename);
      unknown_chromosomes.emplace(chrom);
      continue;
    }
    if (targeted_id == -2) {
      throw error::Error("Could not parse header of bam '{}'", bam_filename);
    }
    if (targeted_id < 0) {
      throw error::Error("Could not determine tid of '{}' in bam '{}'", chrom, bam_filename);
    }
    regions[chrom].emplace_back(start, end);
    ++region_count;
  }
  // Sort the regions for each chromosome
  for (auto& [chromosome, intervals] : regions) {
    if (!std::ranges::is_sorted(intervals)) {
      std::ranges::sort(intervals);
    }
  }
  Logging::Info("Read {} regions from BED file '{}'", region_count, bed_filename.string());
  return regions;
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
