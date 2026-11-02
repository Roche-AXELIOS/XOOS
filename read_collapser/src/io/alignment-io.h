
#pragma once

#include <gsl/gsl>  // NOLINT
#include <optional>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/io/alignment-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "clustering/clustering.h"
#include "core/read-collapser-options.h"
#include "core/region.h"
#include "io/alignment.h"
#include "metrics/metrics.h"

namespace xoos::read_collapser {

/// Returns true if the record has a YF:i:1 BAM auxiliary tag, indicating a rescued secondary alignment
/// from pangenome alignment that should be treated as a primary alignment for clustering.
bool HasRescuedSecondaryTag(const bam1_t* record);

/**
 * Determine if the record in @p alignment_ptr should be discarded based on @p options like mapping quality and
 * alignment. Specifically, the following criteria are applied:
 * 1. If the mapping quality of the record is less than the minimum mapping quality specified
 * 2. If any of the BAM flags specified in options.exclude_flags are set in the record's flags, with a special case for
 *    secondary alignments with YF:i:1 tag when options.markdup_options.cluster_rescued_secondaries is true. In this
 *    case, the read is not discarded if BAM_FSECONDARY is the only matching exclude flag and the YF:i:1 tag is present,
 *    indicating that the read is a rescued secondary alignment that should be clustered.
 * 3. If the discordant duplex error rate is above the threshold specified.
 * 4. If clustering by UMI is enabled and the read is missing both 5' and 3' UMIs.
 * The appropriate metrics in @p metrics are incremented based on the reason for discarding the record.
 */
bool ShouldDiscardRecord(const ReadCollapserOptions& options, const Alignment* alignment_ptr, Metrics& metrics);

/**
 * Read all alignments from @p reader that overlap with @p region, reads will be filtered according to the @p options.
 * The alignments are stored in @p alignments. Discard metrics are accumulated into @p metrics.
 */
void ReadAlignments(const ReadCollapserOptions& options,
                    const io::AlignmentReader& reader,
                    const Region& region,
                    vec<AlignmentPtr>& alignments,
                    Metrics& metrics);

/**
 * HtsIdxWriter will keep the filename of the index in scope to be referenced when the writer is deleted
 * and the index is saved. This is required because htslib will not manage the memory of the filename,
 * and will only keep a non-owned pointer to the filename.
 */
struct HtsIdxWriter {
  /// A non-owning pointer to the BAM file that is being indexed, must remain in scope until after the index is saved.
  htsFile* bam{};
  /// The filename of the index to be saved.
  std::string idx_fn;

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
 * Discard and unmapped-read metrics are accumulated into @p metrics.
 */
void ReadAndWriteUnmappedReads(const ReadCollapserOptions& options,
                               const CreateClusterId& create_cluster_id,
                               const io::AlignmentReader& reader,
                               const UnmappedAlignmentWriter& writer,
                               Metrics& metrics);

/// Remove the PG auxiliary tag from @p record if present.
void RemovePgTag(bam1_t* record);

/**
 * Write a single alignment @p record to @p output. If @p mark_duplicate is true, the alignment is marked as a
 * duplicate. The cluster ID and cluster size are not added to the BAM record.
 */
void WriteAlignment(const bam1_t* record, bool mark_duplicate, const io::HtsFilePtr& output);

/**
 * Write a single alignment @p record to @p output. If @p mark_duplicate is true, the alignment is marked as a
 * duplicate. The cluster ID and cluster size are not added to the BAM record.
 * @param remove_pg_tag If true, strip the PG auxiliary tag from the output record.
 */
void WriteAlignment(const bam1_t* record, bool mark_duplicate, const io::HtsFilePtr& output, bool remove_pg_tag);

/**
 * Write a single alignment @p record to @p output. If @p mark_duplicate is true, the alignment is marked as a
 * duplicate. The cluster ID and cluster size are added to the BAM record.
 */
void WriteAlignment(const bam1_t* record,
                    bool mark_duplicate,
                    const ClusterId& cluster_id,
                    u32 cluster_size,
                    const io::HtsFilePtr& output);

/**
 * Write a single alignment @p record to @p output. If @p mark_duplicate is true, the alignment is marked as a
 * duplicate. The cluster ID and cluster size are added to the BAM record.
 * @param remove_pg_tag If true, strip the PG auxiliary tag from the output record.
 */
void WriteAlignment(const bam1_t* record,
                    bool mark_duplicate,
                    const ClusterId& cluster_id,
                    u32 cluster_size,
                    const io::HtsFilePtr& output,
                    bool remove_pg_tag);

/**
 * Write all alignments in @p alignments to @p output. If @p include_cluster_info is true, the cluster ID and cluster
 * size are added to the BAM records. The @p remove_duplicates determines if duplicates are marked or removed.
 */
void WriteAlignments(const vec<AlignmentPtr>& alignments,
                     bool remove_duplicates,
                     bool include_cluster_info,
                     const io::HtsFilePtr& output);

/**
 * Write all alignments in @p alignments to @p output. If @p include_cluster_info is true, the cluster ID and cluster
 * size are added to the BAM records. The @p remove_duplicates determines if duplicates are marked or removed.
 * @param remove_pg_tag If true, strip the PG auxiliary tag from each output record.
 */
void WriteAlignments(const vec<AlignmentPtr>& alignments,
                     bool remove_duplicates,
                     bool include_cluster_info,
                     const io::HtsFilePtr& output,
                     bool remove_pg_tag);

void ConcatenateBamFiles(const vec<fs::path>& input_bam_files, const fs::path& output_bam_file);

}  // namespace xoos::read_collapser
