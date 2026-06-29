#include "io/alignment-io.h"

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>

#include "bam_cat.h"
#include "util/duplex-util.h"

namespace xoos::read_collapser {

AlignmentReader OpenAlignmentReader(const fs::path& location) {
  auto bam = io::HtsOpen(location, "rb");
  if (bam->format.format != htsExactFormat::sam && bam->format.format != htsExactFormat::bam) {
    throw std::runtime_error("hts_open failed");
  }

  auto idx = io::HtsIdxLoad(location, HTS_FMT_BAI);
  auto hdr = io::SamHdrRead(bam.get());
  return AlignmentReader{std::move(bam), std::move(hdr), std::move(idx)};
}

vec<AlignmentReader> OpenAlignmentReaders(const fs::path& bam_filename, u32 reader_count) {
  vec<AlignmentReader> alignment_readers;
  alignment_readers.reserve(reader_count);
  for (u32 i = 0; i < reader_count; ++i) {
    alignment_readers.emplace_back(OpenAlignmentReader(bam_filename));
  }
  return alignment_readers;
}

// Determine if the alignment in @p record should be discarded based on @p options like mapping quality and alignment
// flags, increment appropriate values of @p metrics, and return true if the alignment should be discarded.
bool ShouldDiscardAlignment(const ReadCollapserOptions& options, const bam1_t* record, Metrics& metrics) {
  bool discard = false;
  ++metrics.input_reads;
  if (record->core.qual < options.min_mapq) {
    ++metrics.discarded_low_mapq_reads;
    discard = true;
  }
  if ((record->core.flag & options.exclude_flags) != 0) {
    ++metrics.discarded_by_flags_reads;
    discard = true;
  }
  if (options.max_discordant_duplex_error_percentage) {
    const auto duplex_error_rate = GetDiscordantDuplexErrorRate(record);
    if (duplex_error_rate && *duplex_error_rate * 100.0 > *options.max_discordant_duplex_error_percentage) {
      ++metrics.discarded_high_discordant_duplex_percentage_reads;
      discard = true;
    }
  }
  if (discard) {
    ++metrics.discarded_total_reads;
  }
  return discard;
}

void ReadAlignments(const ReadCollapserOptions& options,
                    const AlignmentReader& reader,
                    const Region& region,
                    vec<AlignmentPtr>& alignments) {
  const auto itr = io::SamItrQueryI(reader.idx.get(), region.tid, region.start, region.end);
  auto& metrics = ConcurrentMetrics::Get();
  while (true) {
    auto record = io::Bam1Ptr{bam_init1()};
    const auto has_more_reads = io::SamItrNext(reader.bam.get(), itr.get(), record.get());
    if (!has_more_reads) {
      break;
    }
    // Skip alignments with a start position outside of the region. This is to avoid double counting
    // reads which start in one region but overlap another region.
    // TODO: This should be improved to handle alignments which start outside of the region,
    //  but overlap the region. As it is now these alignments will not be included. This will have a
    //  very small affect on some target panels, where some alignments that start outside of the region
    //  but within the wiggle room are not included.
    if (record->core.pos < region.start || record->core.pos >= region.end) {
      continue;
    }
    if (ShouldDiscardAlignment(options, record.get(), metrics)) {
      continue;
    }
    alignments.emplace_back(std::make_shared<Alignment>(std::move(record)));
  }
}

HtsIdxWriter::HtsIdxWriter(HtsIdxWriter&& other) noexcept : bam(other.bam), idx_fn(std::move(other.idx_fn)) {
  // Nullify the other bam pointer to prevent double deletion.
  other.bam = nullptr;
}

HtsIdxWriter& HtsIdxWriter::operator=(HtsIdxWriter&& other) noexcept {
  if (this != &other) {
    // First, clean up any existing resources.
    if (bam != nullptr) {
      io::HtsFlush(bam);
      io::SamIdxSave(bam);
    }
    // Move the resources from other to this.
    bam = other.bam;
    idx_fn = std::move(other.idx_fn);
    // Nullify the other bam pointer to prevent double deletion.
    other.bam = nullptr;
  }
  return *this;
}

HtsIdxWriter::~HtsIdxWriter() {
  // Add a flush first to update the timestamp of the BAM if necessary,
  // this helps to avoid a samtools warning about the index being older than the BAM.
  if (bam != nullptr) {
    io::HtsFlush(bam);
    io::SamIdxSave(bam);
  }
}

AlignmentWriter OpenAlignmentWriter(const fs::path& location, sam_hdr_t* hdr, const bool write_index) {
  htsFormat fmt;
  fmt.category = sequence_data;
  fmt.format = bam;
  fmt.compression = bgzf;
  fmt.compression_level = 1;
  fmt.specific = nullptr;

  create_directories(location.parent_path());
  auto bam = io::HtsOpenFormat(location, "wbz1", &fmt);
  std::optional<HtsIdxWriter> bam_idx;

  if (write_index) {
    bam_idx.emplace(bam.get(), location.string() + ".bai");
    io::SamIdxInit(bam.get(), hdr, 0, bam_idx.value().idx_fn);
  }

  io::SamHdrWrite(bam.get(), hdr);

  return AlignmentWriter{std::move(bam), std::move(bam_idx)};
}

void WriteAlignment(const bam1_t* record, bool mark_duplicate, const io::HtsFilePtr& output) {
  auto record_modified = io::Bam1Ptr{bam_dup1(record)};
  if (mark_duplicate) {
    record_modified->core.flag |= BAM_FDUP;
  }
  io::SamWrite1(output.get(), output->bam_header, record_modified.get());
}

void WriteAlignment(const bam1_t* record,
                    bool mark_duplicate,
                    const ClusterId& cluster_id,
                    u32 cluster_size,
                    const io::HtsFilePtr& output) {
  auto record_modified = io::Bam1Ptr{bam_dup1(record)};
  if (mark_duplicate) {
    record_modified->core.flag |= BAM_FDUP;
  }
  io::BamAuxAppend(record_modified.get(), "DI", fmt::format("{}:{}", cluster_id.super_region_id, cluster_id.count));
  io::BamAuxAppend(record_modified.get(), "DS", cluster_size);
  io::SamWrite1(output.get(), output->bam_header, record_modified.get());
}

void WriteAlignments(const vec<AlignmentPtr>& alignments,
                     const bool remove_duplicates,
                     const bool include_cluster_info,
                     const io::HtsFilePtr& output) {
  for (const auto& alignment : alignments) {
    if (remove_duplicates && alignment->is_duplicate) {
      continue;
    }
    const auto mark_duplicate = !remove_duplicates && alignment->is_duplicate;

    if (auto cluster = alignment->cluster; include_cluster_info && cluster != nullptr) {
      auto cluster_size = static_cast<u32>(cluster->alignments.size());
      WriteAlignment(alignment->record.get(), mark_duplicate, cluster->cluster_id, cluster_size, output);
    } else {
      WriteAlignment(alignment->record.get(), mark_duplicate, output);
    }
  }
}

void ReadAndWriteUnmappedAlignments(const ReadCollapserOptions& options,
                                    const CreateClusterId& create_cluster_id,
                                    const AlignmentReader& reader,
                                    const UnmappedAlignmentWriter& writer) {
  const auto itr = io::HtsItrPtr{sam_itr_queryi(reader.idx.get(), HTS_IDX_NOCOOR, 0, 0)};
  if (itr == nullptr) {
    // If the region is unmapped and the bam file does not contain any reads, itr will be nullptr.
    // Otherwise, as long as the bam file contains at least one read (mapped or unmapped), itr will be a valid pointer
    // for an unmapped region. This is a limitation of sam_itr_queryi, which has special handling of HTS_IDX_NOCOOR and
    // assumes all unmapped reads can be found at the end of the bam file (i.e. when a bam file is empty, this section
    // cannot exist).
    return;
  }
  auto& metrics = ConcurrentMetrics::Get();
  auto record = io::Bam1Ptr{bam_init1()};
  while (io::SamItrNext(reader.bam.get(), itr.get(), record.get())) {
    if (ShouldDiscardAlignment(options, record.get(), metrics)) {
      continue;
    }
    ++metrics.unmapped_reads;
    writer(create_cluster_id(), record.get());
  }
}

void ConcatenateBamFiles(const vec<fs::path>& input_bam_files, const fs::path& output_bam_file) {
  const char* first_file = input_bam_files[0].c_str();
  vec<const char*> input_filenames;
  for (auto& input_bam_file : input_bam_files) {
    input_filenames.emplace_back(input_bam_file.c_str());
  }
  const auto file_count = static_cast<s32>(input_bam_files.size());
  // bam_cat will take ownership of first_sam_file and header
  samFile* first_sam_file = io::HtsOpen(first_file, "r").release();
  const char* output_file = output_bam_file.c_str();
  const s32 bam_cat_result =
      bam_cat(first_sam_file, file_count, input_filenames.data(), nullptr, output_file, nullptr, 1);
  if (bam_cat_result < 0) {
    throw error::Error("bam_cat failed");
  }
}

}  // namespace xoos::read_collapser
