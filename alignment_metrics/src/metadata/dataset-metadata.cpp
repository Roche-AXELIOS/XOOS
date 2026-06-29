#include "metadata/dataset-metadata.h"

#include <filesystem>
#include <locale>
#include <string>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/types/int.h>
#include <xoos/util/string-functions.h>

namespace xoos::alignment_metrics {

// This function checks if a string contains only digits
static bool OnlyDigits(const std::string& str) {
  return std::ranges::all_of(str, [](const char c) { return std::isdigit(c, std::locale::classic()); });
}

/**
 * @brief Extracts dataset-level metadata from a single BAM alignment record.
 *
 * Analyzes the read name format and BAM tags to determine the characteristics of the dataset,
 * including whether it contains post-consensus reads with cluster information, pre-consensus
 * reads with strand/UMI information, and whether it's a duplex or simplex dataset. The function
 * inspects two read name formats:
 *
 * Post-consensus: {cluster_id}-{fwd_partial}-{rev_partial}-{fwd_full}-{rev_full}-{avg_family_size}
 * Pre-consensus: {read_name}|{sid}|{umi_5p}|{umi_3p} or {read_name}|{sid_5p}|{sid_3p}|{umi_5p}|{umi_3p}
 *
 * The presence of the YC tag indicates a duplex dataset. This metadata is used to determine which
 * metrics can be calculated and how to interpret the alignment data throughout the pipeline.
 *
 * @param alignment Smart pointer to a BAM alignment record to analyze
 * @return DatasetMetadata containing flags indicating available dataset features
 */
DatasetMetadata GetDatasetMetadataFromAlignment(const io::Bam1Ptr& alignment) {
  DatasetMetadata metadata;
  const std::string read_name = bam_get_qname(alignment.get());
  auto read_name_parts = string::Split(read_name, "-");
  if (read_name_parts.size() == 6) {
    // If the read name is of the format
    // `{cluster_id}-{forward_partial_count}-{reverse_partial_count}-{forward_full_count}-{reverse_full_count}-{average_family_size}`
    // then it is a post-intermolecular-consensus read, and we
    // can extract the cluster information
    bool last_five_parts_are_numbers = true;
    for (size_t i = 1; i <= 5; ++i) {
      if (!OnlyDigits(read_name_parts[read_name_parts.size() - i])) {
        last_five_parts_are_numbers = false;
        break;
      }
    }
    if (last_five_parts_are_numbers) {
      metadata.has_cluster_info = true;
      // Since it is post-intermolecular-consensus read, there is
      // no more individual strand information
      metadata.has_strand_info = false;
    }
  }
  read_name_parts = string::Split(read_name, "|");
  if (
      // If the read name is of the format for legacy read name format
      // `{read_name}|{sid}|{umi_5p}|{umi_3p}`
      // or
      // `{read_name}|{sid_5p}|{sid_3p}|{umi_5p}|{umi_3p}`
      read_name_parts.size() >= 3 ||
      // Newer read names will use this format, which has the same number of parts when split the post-trimming sid/umi
      // with '|' and the sid from the umi with ':' if the umi is populated
      // `{read_name}|{sid_id}:{umi_5p}:{umi_3p}`
      (read_name_parts.size() == 2 && string::Split(read_name_parts[1], ":").size() == 3)) {
    // umi is found, so it is a pre-intermolecular-consensus read, and we
    // can extract individual strand information as well as
    // read type (full or partial) information
    metadata.has_strand_info = true;
    metadata.has_read_type_info = true;
  }

  // NOTE: this may go away in the future as datasets will have a mixture of reads with and without YC tags
  // Determine if the dataset is duplex or simplex based on the YC tag
  const u8* yc_tag = bam_aux_get(alignment.get(), "YC");
  // the yc_tag is present if the dataset is duplex, if not then it is a simplex dataset
  metadata.is_duplex_dataset = (yc_tag != nullptr);  // NOLINT
  return metadata;
}

/**
 * @brief Determines dataset metadata by reading the first valid alignment from a BAM file.
 *
 * Opens the BAM file, reads the first available alignment record, and extracts dataset metadata
 * from it. Returns a default metadata object if no alignments are found in the file.
 *
 * @param location Path to the BAM file to analyze
 * @return DatasetMetadata extracted from the first alignment, or default values if file is empty
 */
DatasetMetadata GetDatasetMetadataFromBam(const std::filesystem::path& location) {
  const auto bam = io::HtsOpen(location, "rb");
  const auto hdr = io::SamHdrRead(bam.get());
  const auto read = io::Bam1Ptr{bam_init1()};
  if (sam_read1(bam.get(), hdr.get(), read.get()) > 0) {
    return GetDatasetMetadataFromAlignment(read);
  }
  // If we reach here, it means we didn't find any alignments
  // in the BAM file, so we return a default metadata object
  return DatasetMetadata{};
}

}  // namespace xoos::alignment_metrics
