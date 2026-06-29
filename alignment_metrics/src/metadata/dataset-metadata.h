#pragma once

#include <filesystem>

#include <xoos/io/htslib-util/htslib-ptr.h>

namespace xoos::alignment_metrics {

/**
 * Metadata describing what information can be inferred from the read name
 * of alignments in a BAM dataset. This is used to determine what metrics
 * to output as certain metrics are only applicable to certain types of datasets.
 *
 * For example, errors_by_cluster_size is only applicable to post-intermolecular-consensus
 * datasets with cluster information while errors_by_read_type is only applicable if
 * the dataset has read type (full/partial) and/or strand information (forward/reverse).
 */
struct DatasetMetadata {
  // Whether the dataset has cluster information including cluster size,
  // cluster strand (if reads contained in the cluster are all forward, reverse, or mixed),
  // and cluster type (if reads contained in the cluster are all full, partial, or mixed).
  bool has_cluster_info{false};
  // Whether the dataset has individual strand information
  // for each read
  bool has_strand_info{true};
  // Whether the dataset has full/partial read type information
  // for each individual read
  bool has_read_type_info{false};
  // Whether the dataset is a duplex (SBX-D or SBX-FAST) dataset
  // This information is currently inferred from the presence of the YC tag
  // but may change in the future
  bool is_duplex_dataset{false};

  bool operator==(const DatasetMetadata& other) const = default;
};

/**
 * Infer dataset metadata from the read name of a single alignment.
 */
DatasetMetadata GetDatasetMetadataFromAlignment(const io::Bam1Ptr& alignment);

/**
 * Infer dataset metadata from a BAM file.
 */
DatasetMetadata GetDatasetMetadataFromBam(const std::filesystem::path& location);

}  // namespace xoos::alignment_metrics
