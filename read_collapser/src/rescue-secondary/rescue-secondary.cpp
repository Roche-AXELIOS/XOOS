#include "rescue-secondary/rescue-secondary.h"

#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>

#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/io/bam/read-group.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/htslib-util/kstring.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "rescue-secondary/rescue-secondary-metrics.h"
#include "rescue-secondary/rescue-secondary-options.h"

namespace xoos::read_collapser::rescue_secondary {

constexpr std::string_view kDefaultOutputBamFilename = "rescued_secondary.bam";
constexpr std::string_view kRescueTag = "YF";
constexpr u32 kRescueTagValue = 1;

// Pre-reserve size to avoid repeated rehashing on large WGS datasets, which causes
// temporary memory spikes that can trigger OOM errors.
constexpr size_t kEstimatedUnmappedPrimaries = 20'000'000;

constexpr s32 kHtsCacheSizeBytes = 8 * 1024 * 1024;  // 8 MiB

// --- BAM flag predicates ---

namespace {

bool IsPrimaryUnmapped(const bam1_t* const record) {
  const u16 flags = record->core.flag;
  return (flags & BAM_FUNMAP) != 0 && (flags & BAM_FSECONDARY) == 0 && (flags & BAM_FSUPPLEMENTARY) == 0;
}

bool IsSecondary(const bam1_t* const record) {
  return (record->core.flag & BAM_FSECONDARY) != 0;
}

bool IsMapped(const bam1_t* const record) {
  return (record->core.flag & BAM_FUNMAP) == 0;
}

bool IsSupplementary(const bam1_t* const record) {
  return (record->core.flag & BAM_FSUPPLEMENTARY) != 0;
}

bool IsPrimaryAlignment(const bam1_t* const record) {
  const u16 flags = record->core.flag;
  return (flags & BAM_FSECONDARY) == 0 && (flags & BAM_FSUPPLEMENTARY) == 0;
}

// Returns true if the alignment score (AS tag) exceeds min_as_ratio * read_length.
// Uses s64 to preserve sign — AS scores can be negative for some aligners.
bool IsHighQualityAlignment(const bam1_t* const record, const f64 min_as_ratio) {
  const auto read_length = static_cast<f64>(record->core.l_qseq);
  if (read_length <= 0) {
    return false;
  }
  const auto as_score = io::BamAuxGet<s64>(record, "AS");
  if (!as_score.has_value()) {
    return false;
  }
  return static_cast<f64>(as_score.value()) / read_length >= min_as_ratio;
}

// --- I/O helpers ---

io::SamFilePtr OpenInput(const std::string& path, const size_t threads) {
  auto in = io::SamOpen(path, "rb");
  if (threads > 1 && (hts_set_threads(in.get(), static_cast<s32>(threads)) < 0)) {
    throw error::Error("Failed to set {} threads on input BAM '{}'", threads, path);
  }
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg,hicpp-vararg) - hts_set_opt is a C variadic htslib API
  hts_set_opt(in.get(), HTS_OPT_CACHE_SIZE, kHtsCacheSizeBytes);
  return in;
}

io::SamFilePtr OpenOutput(const fs::path& path, const size_t threads) {
  auto out = io::SamOpen(path, "wb");
  if (threads > 1 && (hts_set_threads(out.get(), static_cast<s32>(threads)) < 0)) {
    throw error::Error("Failed to set {} threads on output BAM '{}'", threads, path);
  }
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg,hicpp-vararg) - hts_set_opt is a C variadic htslib API
  hts_set_opt(out.get(), HTS_OPT_CACHE_SIZE, kHtsCacheSizeBytes);
  return out;
}

void AddPgHeader(sam_hdr_t* hdr, const RescueSecondaryOptions& options) {
  const auto& info = options.command_line_info;
  io::SamHdrAddPgLine(
      hdr,
      {.id = info.name, .program_name = info.name, .version_number = info.version, .command_line = info.command_line});
}

bool HasProgramHeader(sam_hdr_t* const hdr, const std::string_view program_name) {
  s32 pos = 0;
  while (const char* const id = sam_hdr_line_name(hdr, io::kPgTag, pos)) {
    if (std::string_view{id} == program_name) {
      return true;
    }

    io::Kstring pn;
    if (sam_hdr_find_tag_id(hdr, io::kPgTag, io::kPgProgramIdTag, id, io::kPgProgramNameTag, pn.Get()) == 0 &&
        pn.Str() == program_name) {
      return true;
    }
    ++pos;
  }
  return false;
}

void ValidateInputNotAlreadyRescued(sam_hdr_t* const hdr, const RescueSecondaryOptions& options) {
  if (HasProgramHeader(hdr, options.command_line_info.name)) {
    throw error::Error(
        "Input BAM was already processed by rescue_secondary; rescue_secondary cannot be run again on it");
  }
}

// Return the first @RG read group from the header, or std::nullopt when none is present.
std::optional<io::ReadGroup> GetFirstReadGroup(sam_hdr_t* const hdr) {
  auto read_groups = io::SamHdrParseReadGroups(hdr);
  if (read_groups.empty()) {
    return std::nullopt;
  }
  return read_groups.front();
}

fs::path BuildOutputBamPath(const RescueSecondaryOptions& options) {
  return file::PrefixedOutputPath(options.output_dir, options.output_prefix, std::string{kDefaultOutputBamFilename});
}

void ValidateOutputDoesNotExist(const fs::path& path, const bool overwrite) {
  if (!overwrite && file::FileExists(path)) {
    throw error::Error("Output file already exists: {}. Use --overwrite to replace.", path);
  }
}

// --- Collated mode processing ---

// Scans buffered secondaries for a query group whose primary was unmapped and
// writes the first high-quality mapped secondary to the output BAM, tagging it
// with YF:i:1. At most one secondary is rescued per query group. Secondaries
// that are unmapped or below the AS/read_length threshold are skipped.
void FlushSecondaries(vec<io::Bam1Ptr>& secondaries,
                      samFile* out,
                      sam_hdr_t* header,
                      RescueSecondaryMetrics& metrics,
                      const f64 min_as_ratio) {
  for (auto& bp : secondaries) {
    bam1_t* rec = bp.get();
    if (IsMapped(rec) && IsHighQualityAlignment(rec, min_as_ratio)) {
      io::BamAuxAppend(rec, std::string{kRescueTag}, kRescueTagValue);
      io::SamWrite1(out, header, rec);
      ++metrics.num_rescued_secondary_records;
      break;  // rescue at most one secondary per query group
    }
  }
}

// Processes a collated BAM stream using inline flushing.
//
// Records are streamed in QNAME order. For each query group:
//   1. The primary alignment is written to output immediately (no copy needed).
//   2. If the primary is unmapped, secondary alignments are buffered via bam_copy1.
//   3. When the next query group starts (or EOF), FlushSecondaries is called to
//      rescue at most one high-quality secondary from the buffer.
//   4. If the primary is mapped, secondaries are only counted — never copied or
//      buffered — which avoids the bam_copy1 overhead for the vast majority of reads.
//   5. Supplementary alignments are always written to the output.
//
// Assumes the primary alignment is the first record in each query name group.
// This ordering is produced when vg giraffe is piped directly into this tool.
void ProcessCollatedBam(const RescueSecondaryOptions& options,
                        const fs::path& output_bam_path,
                        RescueSecondaryMetrics& metrics) {
  const auto in = OpenInput(options.bam, options.threads);
  const auto header = io::SamHdrRead(in.get());
  ValidateInputNotAlreadyRescued(header.get(), options);
  metrics.read_group = GetFirstReadGroup(header.get());
  AddPgHeader(header.get(), options);

  const auto out = OpenOutput(output_bam_path, options.threads);
  io::SamHdrWrite(out.get(), header.get());

  io::Bam1Ptr record{bam_init1()};
  if (!record) {
    throw error::Error("Failed to allocate BAM record");
  }

  // Secondaries are only buffered for query groups with an unmapped primary.
  vec<io::Bam1Ptr> buffered_secondaries;
  std::string current_qname;
  bool primary_unmapped = false;
  bool first_alignment = true;

  // Called at each query group boundary and at EOF. If the primary was unmapped,
  // attempts to rescue one secondary from the buffer; then clears the buffer
  // regardless so it is ready for the next group.
  const auto flush_group = [&]() {
    if (primary_unmapped) {
      FlushSecondaries(buffered_secondaries, out.get(), header.get(), metrics, options.min_alignment_score_ratio);
    }
    buffered_secondaries.clear();
  };

  s32 ret = 0;
  while ((ret = sam_read1(in.get(), header.get(), record.get())) >= 0) {
    ++metrics.num_records;
    const char* qname = bam_get_qname(record.get());

    if (first_alignment) {
      current_qname = qname;
      first_alignment = false;
    }

    // New query group — flush the previous one before starting the next
    if (current_qname != qname) {
      flush_group();
      current_qname = qname;
      primary_unmapped = false;
    }

    bam1_t* rec = record.get();

    if (IsPrimaryAlignment(rec)) {
      // Write the primary immediately — no need to buffer it.
      primary_unmapped = IsPrimaryUnmapped(rec);
      io::SamWrite1(out.get(), header.get(), rec);
      ++metrics.num_primary_records;
    } else if (IsSupplementary(rec)) {
      io::SamWrite1(out.get(), header.get(), rec);
      ++metrics.num_supplementary_records;
    } else if (IsSecondary(rec)) {
      ++metrics.num_secondary_records;
      // Only buffer secondaries when the primary is unmapped. For mapped primaries
      // (the common case), secondaries are counted but never copied, avoiding the
      // per-record bam_copy1 overhead.
      if (primary_unmapped) {
        io::Bam1Ptr rec_copy{bam_init1()};
        if (!rec_copy || bam_copy1(rec_copy.get(), rec) == nullptr) {
          throw error::Error("Failed to copy BAM record for read '{}'", qname);
        }
        buffered_secondaries.push_back(std::move(rec_copy));
      }
    }
  }

  // Flush the last query group (no subsequent QNAME change triggers it)
  flush_group();

  if (ret < -1) {
    throw error::Error("Error reading BAM record from '{}'", options.bam);
  }
}

/**
 * @brief Process a non collated BAM to retain primaries and select secondaries
 *        whose primary alignment is unmapped, tagging them with "YF:i:1"
 *
 * This function performs two passes over the input BAM. Pass 1 collects read names
 * of unmapped primary alignments. Pass 2 writes all primary alignments to the output
 * and writes secondary alignments only if 1) their primary is unmapped, 2) it is
 * the first mapped secondary encountered for that read, and 3) normalized alignment score meets the minimum ratio.
 * Such secondaries receive the "YF:i:1" tag. Supplementary alignments are preserved in the output.
 *
 * @param options RescueSecondaryOptions containing input BAM path and other settings.
 * @param output_bam_path Path to the output BAM file.
 * @param metrics RescueSecondaryMetrics object to store metrics.
 * @throws std::runtime_error On input/output/header/allocation/write failures.
 */
void ProcessNonCollatedBam(const RescueSecondaryOptions& options,
                           const fs::path& output_bam_path,
                           RescueSecondaryMetrics& metrics) {
  // Pass 1: collect read names of unmapped primaries
  Logging::Info("Pass 1: Collecting unmapped primary record read names...");

  std::unordered_set<std::string, std::hash<std::string_view>, std::equal_to<>> unmapped_primary_names;
  unmapped_primary_names.reserve(kEstimatedUnmappedPrimaries);

  {
    const auto in = OpenInput(options.bam, options.threads);
    const auto header = io::SamHdrRead(in.get());
    ValidateInputNotAlreadyRescued(header.get(), options);
    io::Bam1Ptr record{bam_init1()};
    if (!record) {
      throw error::Error("Failed to allocate BAM record");
    }

    s32 ret = 0;
    while ((ret = sam_read1(in.get(), header.get(), record.get())) >= 0) {
      if (IsPrimaryUnmapped(record.get())) {
        const char* qname = bam_get_qname(record.get());
        const auto [it, inserted] = unmapped_primary_names.emplace(qname);
        if (!inserted) {
          throw error::Error("Read '{}' has multiple unmapped primary alignments", qname);
        }
      }
    }
    if (ret < -1) {
      throw error::Error("Error reading BAM record from '{}'", options.bam);
    }
  }

  Logging::Info("Found {} unmapped primary records", unmapped_primary_names.size());

  // Pass 2: write primaries and rescue secondaries
  Logging::Info("Pass 2: Writing output...");

  const auto in = OpenInput(options.bam, options.threads);
  const auto header = io::SamHdrRead(in.get());
  metrics.read_group = GetFirstReadGroup(header.get());
  AddPgHeader(header.get(), options);

  const auto out = OpenOutput(output_bam_path, options.threads);
  io::SamHdrWrite(out.get(), header.get());

  io::Bam1Ptr record{bam_init1()};
  if (!record) {
    throw error::Error("Failed to allocate BAM record");
  }

  s32 ret2 = 0;
  while ((ret2 = sam_read1(in.get(), header.get(), record.get())) >= 0) {
    ++metrics.num_records;
    bam1_t* rec = record.get();

    if (IsPrimaryAlignment(rec)) {
      io::SamWrite1(out.get(), header.get(), rec);
      ++metrics.num_primary_records;
      continue;
    }

    if (IsSupplementary(rec)) {
      io::SamWrite1(out.get(), header.get(), rec);
      ++metrics.num_supplementary_records;
      continue;
    }

    if (!IsSecondary(rec)) {
      continue;
    }

    ++metrics.num_secondary_records;

    if (!IsMapped(rec)) {
      continue;
    }

    // Check the cheap hash lookup first — only a fraction of reads have unmapped primaries.
    // This avoids the more expensive AS tag lookup + float division for the majority of records.
    const char* qname = bam_get_qname(rec);
    const auto it = unmapped_primary_names.find(qname);
    if (it != unmapped_primary_names.end() && IsHighQualityAlignment(rec, options.min_alignment_score_ratio)) {
      io::BamAuxAppend(rec, std::string{kRescueTag}, kRescueTagValue);
      io::SamWrite1(out.get(), header.get(), rec);
      ++metrics.num_rescued_secondary_records;
      // Erase so we rescue at most one secondary per read, and avoid a second hash set.
      unmapped_primary_names.erase(it);
    }
  }
  if (ret2 < -1) {
    throw error::Error("Error reading BAM record from '{}'", options.bam);
  }
}

}  // namespace

void RescueSecondary(const RescueSecondaryOptions& options) {
  file::CreateWritableDirectory(options.output_dir);
  const auto output_bam_path = BuildOutputBamPath(options);
  ValidateOutputDoesNotExist(output_bam_path, options.overwrite);

  const auto metrics_dir = options.output_dir / "metrics";
  file::CreateWritableDirectory(metrics_dir);
  const auto metrics_path =
      file::PrefixedOutputPath(metrics_dir, options.output_prefix, std::string{kDefaultMetricsFilename});
  ValidateOutputDoesNotExist(metrics_path, options.overwrite);

  if (options.bam == kStdin) {
    if (!options.collated) {
      throw error::Error("stdin input requires --collated, as two-pass processing is not possible on a stream");
    }
    Logging::Info("Reading from standard input...");
  } else {
    Logging::Info("Reading from input BAM file: {}", options.bam);
  }

  RescueSecondaryMetrics metrics;

  if (options.collated) {
    Logging::Info("Processing in collated mode (single pass)...");
    ProcessCollatedBam(options, output_bam_path, metrics);
  } else {
    Logging::Info("Processing in non-collated mode (two pass)...");
    ProcessNonCollatedBam(options, output_bam_path, metrics);
  }

  Logging::Info("Wrote {} primary, {} supplementary, and {} rescued secondary records",
                metrics.num_primary_records,
                metrics.num_supplementary_records,
                metrics.num_rescued_secondary_records);

  Logging::Info("Writing rescue secondary metrics to: {}", metrics_path.string());
  const io::TsvMetadata tsv_metadata{.command_line_info = options.command_line_info, .read_group = metrics.read_group};
  metrics.WriteToTsv(metrics_path, tsv_metadata);
}

}  // namespace xoos::read_collapser::rescue_secondary
