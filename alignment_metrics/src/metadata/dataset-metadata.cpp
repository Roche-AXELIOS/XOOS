#include "metadata/dataset-metadata.h"

#include <exception>
#include <filesystem>
#include <string_view>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/read-name/consensus-read-name.h>
#include <xoos/read-name/read-name-parser.h>
#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

/**
 * @brief Extracts dataset-level metadata from a single BAM alignment record.
 *
 * Analyzes the read name format (via the shared `xoos::read_name` parser) and BAM tags to determine
 * the characteristics of the dataset, including whether it contains post-consensus reads with
 * cluster information, pre-consensus reads with strand information, and whether individual reads can
 * be classified as full/partial. The function recognizes:
 *
 * Post-consensus: {cluster_id}-{fwd_partial}-{rev_partial}-{fwd_full}-{rev_full}-{avg_family_size}
 * Pre-consensus: the shared SBX read-name format (legacy `<info>|<sid>|<umi5p>|<umi3p>` or newer
 *   `<info>:<bitflag>|<sid>[:<umi5p>:<umi3p>]`).
 *
 * Read-type (full/partial) information is available whenever the read name carries decodable
 * adapter/UMI information — i.e. the SID/UMI presence bitfield of the newer format or the UMI fields
 * of the legacy format. This means non-UMI datasets that encode the adapter bitfield also expose
 * read-type information, not just UMI datasets.
 *
 * The presence of the YC tag indicates a duplex dataset. This metadata is used to determine which
 * metrics can be calculated and how to interpret the alignment data throughout the pipeline.
 *
 * @param alignment Smart pointer to a BAM alignment record to analyze
 * @return DatasetMetadata containing flags indicating available dataset features
 */
DatasetMetadata GetDatasetMetadataFromAlignment(const io::Bam1Ptr& alignment) {
  DatasetMetadata metadata;
  const std::string_view read_name = bam_get_qname(alignment.get());
  if (read_name::IsConsensusReadName(read_name)) {
    // Post-intermolecular-consensus read: cluster information is available and there is no more
    // individual per-read strand information.
    metadata.has_cluster_info = true;
    metadata.has_strand_info = false;
  } else {
    // For pre-consensus reads, strand comes from the BAM flag and is always available. Read-type
    // (full/partial) information is available whenever the read name carries decodable adapter/UMI
    // information — including non-UMI reads that encode the SID presence bitfield.
    try {
      const bool expect_umis = read_name::HasUmiFields(read_name);
      // Parsing succeeds only when the name carries adapter/UMI information we can decode.
      read_name::ParseReadName(read_name, expect_umis);
      metadata.has_read_type_info = true;
    } catch (const std::exception&) {
      metadata.has_read_type_info = false;
    }
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
 * from it. Returns a default metadata object with `is_empty` set if no alignments are found in the file.
 *
 * @param location Path to the BAM file to analyze
 * @return DatasetMetadata extracted from the first alignment, or default values with `is_empty` set if the
 *   file has zero alignment records
 */
DatasetMetadata GetDatasetMetadataFromBam(const std::filesystem::path& location) {
  const auto bam = io::HtsOpen(location, "rb");
  const auto hdr = io::SamHdrRead(bam.get());
  const auto read = io::Bam1Ptr{bam_init1()};
  if (io::SamRead1(bam.get(), hdr.get(), read.get())) {
    return GetDatasetMetadataFromAlignment(read);
  }
  // If we reach here, it means we didn't find any alignments
  // in the BAM file, so we return a default metadata object flagged as empty
  DatasetMetadata metadata{};
  metadata.is_empty = true;
  return metadata;
}

}  // namespace xoos::alignment_metrics
