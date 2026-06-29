
#pragma once

#include <gsl/gsl>  // NOLINT

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "clustering/clustering.h"
#include "core/read-collapser-options.h"
#include "core/region.h"
#include "io/alignment.h"
#include "metrics/metrics.h"

namespace xoos::read_collapser {

struct AlignmentReader {
  io::HtsFilePtr bam;
  io::SamHdrPtr hdr;
  io::HtsIdxPtr idx;
};

AlignmentReader OpenAlignmentReader(const fs::path& location);

vec<AlignmentReader> OpenAlignmentReaders(const fs::path& bam_filename, u32 reader_count);

bool ShouldDiscardAlignment(const ReadCollapserOptions& options, const bam1_t* record, Metrics& metrics);

/**
 * Read all alignments from @p reader that overlap with @p region, reads will be filtered according to the @p options.
 * The alignments are stored in @p alignments.
 */
void ReadAlignments(const ReadCollapserOptions& options,
                    const AlignmentReader& reader,
                    const Region& region,
                    vec<AlignmentPtr>& alignments);

/**
 * HtsIdxWriter will keep the filename of the index in scope to be referenced when the writer is deleted
 * and the index is saved. This is required because htslib will not manage the memory of the filename,
 * and will only keep a non-owned pointer to the filename.
 */
struct HtsIdxWriter {
  /// A non-owning pointer to the BAM file that is being indexed, must remain in scope until after the index is saved.
  htsFile* bam{};
  /// The filename of the index to be saved.
  std::string idx_fn{};

  HtsIdxWriter(htsFile* const bam_ptr, std::string filename) noexcept : bam(bam_ptr), idx_fn(std::move(filename)) {
  }

  HtsIdxWriter(HtsIdxWriter&& other) noexcept;
  HtsIdxWriter(const HtsIdxWriter&) = delete;

  HtsIdxWriter& operator=(HtsIdxWriter&& other) noexcept;
  HtsIdxWriter& operator=(const HtsIdxWriter&) = delete;

  // Custom destructor to flush bam and save index.
  ~HtsIdxWriter();
};

struct AlignmentWriter {
  // The order of members is important here, bam must be destroyed before bam_idx
  io::HtsFilePtr bam;
  std::optional<HtsIdxWriter> bam_idx;
};

AlignmentWriter OpenAlignmentWriter(const fs::path& location, sam_hdr_t* hdr, bool write_index = true);

using UnmappedAlignmentWriter = std::function<void(const ClusterId&, const bam1_t*)>;

/**
 * Read all unmapped alignments from @p reader and write them to @p writer. The alignments are written with a cluster
 * size of 1. This is used to handle unmapped reads in the BAM file. The cluster ID is generated with @p
 * create_cluster_id for each unmapped read to ensure that each unmapped read has a unique cluster ID.
 */
void ReadAndWriteUnmappedAlignments(const ReadCollapserOptions& options,
                                    const CreateClusterId& create_cluster_id,
                                    const AlignmentReader& reader,
                                    const UnmappedAlignmentWriter& writer);

/**
 * Write a single alignment @p record to @p output. If @p mark_duplicate is true, the alignment is marked as a
 * duplicate. The cluster ID and cluster size are not added to the BAM record.
 */
void WriteAlignment(const bam1_t* record, bool mark_duplicate, const io::HtsFilePtr& output);

/**
 * Write a single alignment @p record to @p output. If @p mark_duplicate is true, the alignment is marked as a
 * duplicate. The cluster ID and cluster size are added to the BAM record if @p cluster_id is not empty.
 */
void WriteAlignment(const bam1_t* record,
                    bool mark_duplicate,
                    const ClusterId& cluster_id,
                    u32 cluster_size,
                    const io::HtsFilePtr& output);
/**
 * Write all alignments in @p alignments to @p output. If @p include_cluster_info is true, the cluster ID and cluster
 * size are added to the BAM records. The @p remove_duplicates determines if duplicates are marked or removed.
 */
void WriteAlignments(const vec<AlignmentPtr>& alignments,
                     bool remove_duplicates,
                     bool include_cluster_info,
                     const io::HtsFilePtr& output);

void ConcatenateBamFiles(const vec<fs::path>& input_bam_files, const fs::path& output_bam_file);

}  // namespace xoos::read_collapser
