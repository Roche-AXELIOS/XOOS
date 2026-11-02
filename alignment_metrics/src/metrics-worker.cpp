#include "metrics-worker.h"

#include <algorithm>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include <ankerl/unordered_dense.h>
#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/io/fasta-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/htslib-util/read-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/util/sequence-functions.h>
#include <xoos/yc-decode/yc-decoder.h>

#include "alignment-metrics-options.h"
#include "core/alignment.h"
#include "core/interval.h"
#include "core/pileup.h"
#include "core/super-region.h"
#include "metrics/accuracy-metrics/hp-stat.h"
#include "metrics/metrics.h"
#include "metrics/stratified-metrics.h"
#include "util/homopolymer-util.h"

namespace xoos::alignment_metrics {

namespace {

void InitializeHomopolymerStats(const std::vector<ReferenceHomopolymer>& hp_regions,
                                const SuperRegion& super_region,
                                StratifiedMetrics& metrics,
                                std::vector<Pileup>& pileups) {
  for (const auto& [start, end, base, labels] : hp_regions) {
    const auto hp_length = static_cast<u8>(end - start + 1);
    const char hp_base = base;
    // Each reference homopolymer carries the label bit vector of its subregion (see
    // FindReferenceHomopolymers). Increment the total HP count only in the label containers this
    // homopolymer belongs to, so per-label HP denominators stay correct.
    metrics.ForEachContainer(labels, [hp_base, hp_length](Metrics& target) {
      if (target.accuracy_metrics.has_value() && target.accuracy_metrics->hp_stats.has_value()) {
        target.accuracy_metrics->hp_stats->IncrementHpCount(HpStatKey{.base = hp_base, .hp_length = hp_length});
      }
    });
    for (auto i = start; i <= end; ++i) {
      if (pileups.at(i - ToUnsigned(super_region.start)).valid) {
        pileups.at(i - ToUnsigned(super_region.start)).is_part_of_hp = true;
        if (i == start) {
          pileups.at(i - ToUnsigned(super_region.start)).reference_context.hp_length = hp_length;
          pileups.at(i - ToUnsigned(super_region.start)).reference_context.hp_base = hp_base;
        }
        pileups.at(i - ToUnsigned(super_region.start)).reference_context.relative_position_within_hp =
            static_cast<u8>(i - start);
      }
    }
  }
}

void InitializeGcWindows(const std::string& ref_seq,
                         const SuperRegion& super_region,
                         const u32 gc_window_size,
                         StratifiedMetrics& metrics,
                         vec<Pileup>& pileups,
                         const bool treat_non_acgt_as_n) {
  constexpr u8 kMaxNCountInGcWindow = 4;
  const auto ref_length = ref_seq.size();
  if (ref_length < gc_window_size || pileups.empty()) {
    return;
  }
  const auto last_window_start = ref_length - gc_window_size;
  // A running count of G/C bases in the current sliding window
  u32 gc_count = 0;
  // A running count of N bases in the current sliding window
  u32 n_count = 0;
  // Note: No need to convert to upper case because it's assumed that the reference sequence
  // is already in upper case (handled by FastaReader).
  const auto is_gc = [](char b) { return b == 'G' || b == 'C'; };
  const auto is_n = [treat_non_acgt_as_n](char b) { return b == 'N' || (treat_non_acgt_as_n && !sequence::IsACGT(b)); };
  // Initialize the first window [0, gc_window_size)
  for (u32 j = 0; j < gc_window_size; ++j) {
    const char base = ref_seq[j];
    gc_count += static_cast<u32>(is_gc(base));
    n_count += static_cast<u32>(is_n(base));
  }
  size_t subregion_index = 0;
  for (size_t i = 0; i <= std::min(last_window_start, pileups.size() - 1); ++i) {
    // Move to the next subregion if the current position exceeds its end
    while (subregion_index < super_region.subregions.size() &&
           ToSigned(i) + super_region.start >= super_region.subregions[subregion_index].end) {
      ++subregion_index;
    }
    if (pileups.at(i).valid && n_count <= kMaxNCountInGcWindow) {
      const auto gc_percentage = static_cast<u8>(std::round((static_cast<f64>(gc_count) / gc_window_size) * 100.0));
      pileups.at(i).reference_context.window_gc_percentage = ToSigned(gc_percentage);
      metrics.ForEachContainer(super_region.subregions.at(subregion_index).labels, [gc_percentage](Metrics& target) {
        if (target.gc_bias_metrics.has_value() && gc_percentage < target.gc_bias_metrics->windows_per_gc_pct.size()) {
          auto& window_count_at_gc_percentage = target.gc_bias_metrics->windows_per_gc_pct.at(gc_percentage);
          ++window_count_at_gc_percentage;
          ++target.gc_bias_metrics->total_windows;
        }
      });
    }
    if (i < last_window_start) {
      // Slide the window: remove the base leaving the window and add the base entering the window
      const char base_leaving = ref_seq[i];
      const char base_entering = ref_seq[i + gc_window_size];
      gc_count -= static_cast<u32>(is_gc(base_leaving));
      n_count -= static_cast<u32>(is_n(base_leaving));
      gc_count += static_cast<u32>(is_gc(base_entering));
      n_count += static_cast<u32>(is_n(base_entering));
    }
  }
}

// Determine whether an alignment should still be processed after trimming. Returns false when trimming has left the
// alignment with no reference span / zero length, or moved it out of every region of interest.
bool ShouldProcessAfterTrim(const Alignment& alignment,
                            const AlignmentMetricsOptions& options,
                            const RegionLookupTable& region_lookup_table,
                            const size_t super_region_index,
                            const size_t subregion_index) {
  if (options.trim_trailing_bases == 0 && options.trim_leading_bases == 0) {
    return true;
  }
  // If the alignment no longer spans any reference positions or its length
  // has become zero after trimming, we skip further processing
  if (alignment.reference_length <= 0 || alignment.read_length <= 0) {
    return false;
  }
  // Trimming may also cause an alignment to fall outside of regions of interest
  // so we need to check if the trimmed alignment is still in any valid region
  return region_lookup_table.IsAlignmentInAnyRegion(
      alignment.rpos, alignment.reference_length, super_region_index, subregion_index);
}

s64 GetGcWindowPositionForRead(const Alignment& alignment, const u32 gc_window_size) {
  const bool is_reverse = (alignment.alignment_ptr->core.flag & BAM_FREVERSE) != 0;
  const s64 reference_end = alignment.rpos + alignment.reference_length;
  const s64 unclamped_anchor = is_reverse ? reference_end - static_cast<s64>(gc_window_size) : alignment.rpos;
  return std::clamp(unclamped_anchor, alignment.rpos, reference_end - 1);
}

std::optional<size_t> GetGcBinIndex(const vec<Pileup>& pileups,
                                    const Alignment& alignment,
                                    const SuperRegion& super_region,
                                    const u32 gc_window_size) {
  const s64 gc_index = GetGcWindowPositionForRead(alignment, gc_window_size) - super_region.start;
  if (gc_index >= 0 && static_cast<size_t>(gc_index) < pileups.size()) {
    const auto& gc_bin = pileups.at(static_cast<size_t>(gc_index)).reference_context.window_gc_percentage;
    if (gc_bin >= 0 && static_cast<size_t>(gc_bin) <= kMaxGcPercentage) {
      return static_cast<size_t>(gc_bin);
    }
  }
  return std::nullopt;
}

bool IsMismatch(const bool include_ambiguous_reference_bases, const char ref_base, const char query_base) {
  if (include_ambiguous_reference_bases && sequence::IsAmbiguousBase(ref_base)) {
    return !sequence::AmbiguousBaseMatchesQuery(ref_base, query_base);
  }
  return query_base != ref_base;
}

bool IsValidReferenceBase(const bool include_ambiguous_reference_bases, const char ref_base) {
  return sequence::IsACGT(ref_base) || (include_ambiguous_reference_bases && sequence::IsAmbiguousBase(ref_base));
}

std::unique_ptr<io::FastaReader> CreateFastaReader(const fs::path& reference_input) {
  if (reference_input.empty()) {
    return nullptr;
  }
  if (!fs::exists(reference_input)) {
    throw error::Error(fmt::format("Reference file does not exist: {}", reference_input.string()));
  }
  return std::make_unique<io::FastaReader>(reference_input.string());
}
}  // namespace

BamItrAlignmentProvider::BamItrAlignmentProvider(const SuperRegion& super_region,
                                                 const io::AlignmentReader& alignment_reader)
    : _alignment_reader(alignment_reader),
      _itr(io::SamItrRegion(
          alignment_reader.idx.get(), alignment_reader.hdr.get(), CreateHtsRegList(super_region).release(), 1)) {
}

s32 BamItrAlignmentProvider::operator()(bam1_t* const bam1_ptr) const {
  return io::SamItrMultiNext(_alignment_reader.bam_file.get(), _itr.get(), bam1_ptr);
}

MetricsWorker::MetricsWorker(const AlignmentMetricsOptions& options,
                             const DatasetMetadata& dataset_metadata,
                             const RegionLookupTable& region_lookup_table)
    : _options(options),
      _dataset_metadata(dataset_metadata),
      _region_lookup_table(region_lookup_table),
      _alignment_reader(!_options.bam_input.empty()
                            ? std::make_shared<io::AlignmentReader>(io::OpenAlignmentReader(_options.bam_input))
                            : nullptr),
      _fasta_reader(CreateFastaReader(options.reference_input)),
      _labeled_metrics(
          std::make_shared<StratifiedMetrics>(_options, _dataset_metadata, _region_lookup_table.GetLabels())),
      // Aliasing shared_ptr sharing ownership with `_labeled_metrics` but pointing at its label-0
      // container, so GetMetrics() and existing single-container tests keep working.
      _metrics(_labeled_metrics, &_labeled_metrics->At(0)) {
}

void MetricsWorker::Initialize(const SuperRegion& super_region,
                               const size_t super_region_id,
                               std::optional<std::string>&& ref_seq) {
  // Reset pileups
  if (std::cmp_less(_pileups.size(), super_region.end - super_region.start)) {
    // If the current pileup size is smaller than the required size, we clear and resize
    // since reallocation may occur anyway and it will not save much time to reset the values
    // in place just to have them reallocated and copied again.
    _pileups.clear();
  } else {
    for (auto& pileup : _pileups) {
      if (pileup.valid) {
        pileup.depth_profile.Reset();
        pileup.error_profile.Reset();
        pileup.hp_error_profile.Reset();
        pileup.reference_context.Reset();
        pileup.valid = false;
        pileup.is_part_of_hp = false;
      }
    }
  }
  _pileups.resize(static_cast<size_t>(super_region.end - super_region.start));
  _ref_seq = std::move(ref_seq);
  // Set valid positions in the pileups
  for (const auto& subregion : super_region.subregions) {
    for (s64 i = subregion.start; i < subregion.end; ++i) {
      const auto pileup_index = ToUnsigned(i - super_region.start);
      const bool has_reference_base = _ref_seq.has_value() && pileup_index < _ref_seq->size();
      _pileups.at(pileup_index).valid =
          !has_reference_base ||
          IsValidReferenceBase(_options.include_ambiguous_reference_bases, _ref_seq->at(pileup_index));
      // Typically, we skip positions that do not have coverage when calculating coverage metrics
      // but if the user has specified to include empty positions, we will mark them as having coverage in the pileup
      // so that they are included in the coverage histogram calculations.
      if (!_options.exclude_empty_positions) {
        _pileups.at(pileup_index).depth_profile.has_coverage = true;
      }
    }
  }
  // Initialize GC windows if GC bias metrics are enabled and a reference sequence is available
  if (_options.metric_types.has_gc_bias_metrics && _ref_seq.has_value()) {
    InitializeGcWindows(*_ref_seq,
                        super_region,
                        _options.gc_scan_window_size,
                        *_labeled_metrics,
                        _pileups,
                        !_options.include_ambiguous_reference_bases);
  }
  // If we are calculating homopolymer metrics, set `is_part_of_hp` to true for all positions in the homopolymer regions
  if (_options.calculate_hp_metrics) {
    const std::vector<ReferenceHomopolymer> hp_regions = FindReferenceHomopolymers(*_ref_seq,
                                                                                   _options.min_hp_length,
                                                                                   _options.max_hp_length,
                                                                                   super_region,
                                                                                   _options.hp_subsampling_fraction,
                                                                                   _options.hp_subsampling_seed);
    InitializeHomopolymerStats(hp_regions, super_region, *_labeled_metrics, _pileups);
  }
  // Clear internal states
  _current_subregion_index = 0;
  _current_super_region = super_region;
  _current_super_region_index = super_region_id;
}

void MetricsWorker::ProcessRegion(const SuperRegion& super_region, const size_t super_region_id) {
  const BamItrAlignmentProvider alignment_provider(super_region, *_alignment_reader);
  // Fetch reference sequence from specified path
  std::optional<std::string> ref_seq;
  if (!_options.reference_input.empty()) {
    try {
      // Fetch the reference sequence for the specified region, ensuring all bases are uppercase.
      ref_seq = _fasta_reader->GetSequence(super_region.chromosome,
                                           super_region.start,
                                           super_region.end + static_cast<s64>(_options.reference_padding),
                                           true);
    } catch (const std::runtime_error& e) {  // NOSONAR cpp:S1181 - GetSequence only throws std::runtime_error
      Logging::Error("Error reading reference sequence for region {}:{}-{}\n{}",
                     super_region.chromosome,
                     super_region.start,
                     super_region.end,
                     e.what());
      // Throw if the reference sequence is missing from the FASTA file or the requested region is not fully covered.
      // If we simply skip the contig, we would miss the reads that are aligned to this contig, leading to incomplete
      // metrics.
      throw;
    }
  }
  ProcessRegion(super_region, super_region_id, alignment_provider, std::move(ref_seq));
}

void MetricsWorker::ProcessNoCoverageRegion(const std::string& chromosome, const Interval& region) const {
  std::optional<std::string> ref_seq;
  if (!_options.reference_input.empty()) {
    try {
      if (_fasta_reader->HasContig(chromosome)) {
        ref_seq = _fasta_reader->GetSequence(
            chromosome, region.start, region.end + static_cast<s64>(_options.reference_padding), true);
      } else {
        // For coverage, we can approximate the size using the BAM header length. The only downside is we cannot
        // exclude the ambiguous (non-ACGT) bases in this contig from the denominator for coverage calculation but the
        // impact on the overall coverage metrics is negligible so we log a warning instead of throwing an error.
        if (_options.metric_types.has_coverage_metrics) {
          Logging::Warn(
              "Contig '{}' is listed in the BAM header but has no mapped reads and is missing from the reference "
              "FASTA. The BAM header length will be used for coverage calculation without accounting for ambiguous "
              "bases.",
              chromosome);
        }
        // For GC bias, we cannot approximate the missing reference sequence because we need the actual sequence for GC
        // content calculation. Therefore, we throw an error if GC bias metrics are enabled and the reference sequence
        // is missing.
        if (_options.metric_types.has_gc_bias_metrics) {
          throw error::Error(
              "Contig '{}' is listed in the BAM header but has no mapped reads and is missing from the reference "
              "FASTA. GC windows cannot be collected from this contig.",
              chromosome);
        }
      }
    } catch (const std::runtime_error& error) {  // NOSONAR cpp:S1181 - GetSequence only throws std::runtime_error
      Logging::Error("Unable to read reference sequence for region {}:{}-{}\n{}",
                     chromosome,
                     region.start,
                     region.end,
                     error.what());
      throw;
    }
  }
  ProcessNoCoverageRegion(region, std::move(ref_seq));
}

void MetricsWorker::ProcessNoCoverageRegion(const Interval& region, std::optional<std::string> ref_seq) const {
  if (!_options.exclude_empty_positions && _options.metric_types.has_coverage_metrics) {
    u64 region_length = 0;
    if (ref_seq.has_value()) {
      const bool include_ambiguous = _options.include_ambiguous_reference_bases;
      // Restrict the region length to only include valid reference bases
      // and make sure to stay within the bounds of both the reference sequence and the region.
      const auto reference_end = ref_seq->begin() + std::min(region.end - region.start, ToSigned(ref_seq->size()));
      region_length = ToUnsigned(std::count_if(ref_seq->begin(), reference_end, [include_ambiguous](const char c) {
        return IsValidReferenceBase(include_ambiguous, c);
      }));
    } else {
      region_length = ToUnsigned(region.end - region.start);
    }
    // Route the empty (zero-coverage) positions into every label container this region belongs to.
    _labeled_metrics->ForEachContainer(region.labels, [region_length](Metrics& target) {
      target.coverage_metrics->AddEmptyHistogramData(region_length);
    });
  }
  // Update the GC window count for the no-coverage region if GC bias metrics are enabled
  if (_options.metric_types.has_gc_bias_metrics && ref_seq.has_value()) {
    // The pileup array is used purely for the purpose of tracking valid/invalid
    // reference positions within the no-coverage region.
    vec<Pileup> pileups(ToUnsigned(region.end - region.start));
    for (size_t i = 0; i < pileups.size(); ++i) {
      if (i < ref_seq->size()) {
        pileups.at(i).valid = IsValidReferenceBase(_options.include_ambiguous_reference_bases, (*ref_seq)[i]);
      }
    }
    // We don't actually make use of the chromosome name here, so we use a placeholder value.
    InitializeGcWindows(
        ref_seq.value(),
        {.chromosome = kUnknownChromosome, .start = region.start, .end = region.end, .subregions = {region}},
        _options.gc_scan_window_size,
        *_labeled_metrics,
        pileups,
        !_options.include_ambiguous_reference_bases);
  }
}

void MetricsWorker::UpdateHomopolymerSpanningCoverage(Alignment& alignment, const s64 ref_pos) {
  const auto ref_end = alignment.rpos + alignment.reference_length;
  for (s64 i = ref_pos; i < alignment.reference_length + alignment.rpos; ++i) {
    const size_t relative_position = i - _current_super_region.start;
    if (relative_position < _pileups.size() && _pileups.at(relative_position).valid &&
        _pileups.at(relative_position).is_part_of_hp &&
        _pileups.at(relative_position).reference_context.relative_position_within_hp == 0) {
      ++_pileups.at(relative_position).hp_error_profile.total_reads;
      // If the read spans the entire HP region with at least `kRequiredAnchorOverlap` bases on each
      // side, then it contributes to the spanning coverage
      const s64 hp_start = i;
      const s64 hp_end = i + _pileups.at(relative_position).reference_context.hp_length;
      if (std::cmp_greater_equal(hp_start - ref_pos, kRequiredAnchorOverlap) &&
          std::cmp_greater_equal(ref_end - hp_end, kRequiredAnchorOverlap)) {
        ++_pileups.at(relative_position).hp_error_profile.spanning_reads;
        // Each read carries an `hp_error_profile` object that records the errors and base qualities
        // for each homopolymer region covered by the read. This object is aggregated into the pileup's
        // `hp_error_profile` object after the read has been fully processed. We initialize the `hp_error_profile`
        // by marking the `valid` field to true for the first base of the homopolymer region that is covered by the
        // read.
        alignment.hp_error_profile[hp_start].valid = true;
      }
    }
  }
}

void MetricsWorker::ProcessAlignment(Alignment& alignment) {
  // Get the reference position and query position for the first CIGAR operation
  // We extract the relevant trimming information from the CIGAR string
  alignment.Trim(_options.trim_leading_bases, _options.trim_trailing_bases);
  // Check if the alignment should still be processed after trimming
  if (!ShouldProcessAfterTrim(
          alignment, _options, _region_lookup_table, _current_super_region_index, _current_subregion_index)) {
    return;
  }
  const bool passes_read_filter = alignment.alignment_ptr->core.qual >= _options.min_mapq &&
                                  (alignment.alignment_ptr->core.flag & _options.exclude_flags) == 0;
  // Check if the alignment is safe to count for read-level metrics. This means that the read
  // does not overlap with other subregions in another super region.
  // If the read is placed unmapped, it does not have a reference length and should only be
  // picked up by one super region (the region containing the rpos) so we consider it safe to count.
  const bool is_safe_to_count =
      _region_lookup_table.IsSafeToCountAlignment(_current_super_region.chromosome,
                                                  {alignment.rpos, alignment.rpos + alignment.reference_length},
                                                  _current_super_region_index,
                                                  _current_subregion_index) ||
      (alignment.alignment_ptr->core.flag & BAM_FUNMAP) != 0;
  // Check if the current alignment is already handled by another subregion in another super region
  // If not we will update read-level metrics for this alignment
  if (is_safe_to_count && _options.metric_types.has_read_metrics) {
    RecordRead(alignment, passes_read_filter);
  }
  if (_options.metric_types.has_gc_bias_metrics) {
    RecordReadForGcBias(alignment, passes_read_filter);
  }

  if (!passes_read_filter) {
    // Skip further processing if the read does not pass the filter
    return;
  }

  s64 ref_pos = alignment.rpos;
  u32 query_pos = alignment.qpos;

  // Increment total coverage and spanning coverage for homopolymer positions that are covered by this read
  if (_options.metric_types.has_accuracy_metrics && _options.calculate_hp_metrics) {
    UpdateHomopolymerSpanningCoverage(alignment, ref_pos);
  }
  // Iterate through each CIGAR operation if base-level metrics (accuracy or coverage) are enabled
  if (_options.metric_types.has_accuracy_metrics || _options.metric_types.has_coverage_metrics) {
    ProcessCigarOps(alignment, ref_pos, query_pos);
  }

  // After processing all CIGAR operations, handle homopolymer errors for the read if applicable
  if (_options.metric_types.has_accuracy_metrics && _options.calculate_hp_metrics) {
    HandleHomopolymerErrors(alignment);
  }
}

void MetricsWorker::RecordRead(const Alignment& alignment, const bool passes_read_filter) const {
  // Determine the labels a read is attributed to via union: it counts toward every label whose
  // region of interest its reference span overlaps. Unmapped reads have no span and no labels; they
  // are routed to all containers below.
  const bool is_unmapped = (alignment.alignment_ptr->core.flag & BAM_FUNMAP) != 0;
  const u16 read_labels = is_unmapped ? 0
                                      : _region_lookup_table.GetLabelsForInterval(
                                            _current_super_region.chromosome,
                                            {alignment.rpos, alignment.rpos + alignment.reference_length});
  // Process read-level metrics before read-level filtering because we count the total number of reads
  // without filtering for some of the read-level metrics.
  const auto add_read = [this, &alignment, passes_read_filter](Metrics& target) {
    if (target.read_metrics.has_value()) {
      target.read_metrics->AddRead(alignment,
                                   alignment.alignment_metadata,
                                   SaTagFilter{.min_mapq = _options.min_mapq, .exclude_flags = _options.exclude_flags},
                                   passes_read_filter,
                                   _region_lookup_table);
    }
  };
  if (is_unmapped) {
    // A placed unmapped read belongs to no labeled region; count it in every label container so each
    // label's read metrics mirror an independent single-label run.
    _labeled_metrics->ForAllContainers(add_read);
  } else {
    _labeled_metrics->ForEachContainer(read_labels, add_read);
  }
}

void MetricsWorker::RecordReadForGcBias(const Alignment& alignment, const bool passes_read_filter) const {
  constexpr u16 kGcBiasExcludeFlags = BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY;
  if (alignment.reference_length <= 0 || (alignment.alignment_ptr->core.flag & kGcBiasExcludeFlags) != 0) {
    return;
  }
  const auto gc_bin_index = GetGcBinIndex(_pileups, alignment, _current_super_region, _options.gc_scan_window_size);
  if (!gc_bin_index.has_value()) {
    return;
  }
  const s64 gc_window_pos = GetGcWindowPositionForRead(alignment, _options.gc_scan_window_size);
  const u16 gc_window_labels =
      _region_lookup_table.GetLabelsForInterval(_current_super_region.chromosome, {gc_window_pos, gc_window_pos + 1});
  _labeled_metrics->ForEachContainer(gc_window_labels, [&gc_bin_index, passes_read_filter](Metrics& target) {
    if (target.gc_bias_metrics.has_value()) {
      target.gc_bias_metrics->AddRead(*gc_bin_index, passes_read_filter);
    }
  });
}

void MetricsWorker::ProcessCigarOps(Alignment& alignment, s64 ref_pos, u32 query_pos) {
  for (size_t cigar_index = 0; cigar_index < alignment.cigars.size(); ++cigar_index) {
    auto [cigar_op, cigar_len] = alignment.GetCigarAt(cigar_index);
    switch (cigar_op) {
      case BAM_CMATCH:
      case BAM_CEQUAL:
      case BAM_CDIFF:
        std::tie(ref_pos, query_pos) = HandleMatchCigar(alignment, ref_pos, query_pos, cigar_len);
        break;
      case BAM_CINS: {
        const bool is_last_cigar_op = (cigar_index == alignment.cigars.size() - 1);
        std::tie(ref_pos, query_pos) = HandleInsertionCigar(alignment, ref_pos, query_pos, cigar_len, is_last_cigar_op);
        break;
      }
      case BAM_CDEL:
      case BAM_CREF_SKIP:
        std::tie(ref_pos, query_pos) = HandleDeletionCigar(alignment, ref_pos, query_pos, cigar_len);
        break;
      case BAM_CSOFT_CLIP: {
        query_pos += cigar_len;
        break;
      }
      case BAM_CHARD_CLIP:
      case BAM_CPAD:
        break;
      default:
        Logging::Warn(
            "Unknown CIGAR operation {} in alignment for {}", cigar_op, bam_get_qname(alignment.alignment_ptr));
        break;
    }
  }
}

void MetricsWorker::RecordMatch(
    Pileup& pileup, const Alignment& alignment, const s64 rpos, const u32 query_pos, const u32 query_offset) {
  if (!_options.metric_types.has_accuracy_metrics && !_options.metric_types.has_coverage_metrics) {
    return;
  }
  const char query_base = io::GetBase(bam_get_seq(alignment.alignment_ptr), query_pos + query_offset);
  if (_options.metric_types.has_accuracy_metrics) {
    // Check if rpos is within the region of interest
    const bool ref_position_valid = _ref_seq.has_value() && !_ref_seq->empty() &&
                                    rpos - _current_super_region.start >= 0 &&
                                    rpos - _current_super_region.start < ToSigned(_ref_seq->size());
    const auto ref_index = ref_position_valid ? ToUnsigned(rpos - _current_super_region.start) : 0;
    const char ref_base = ref_position_valid ? _ref_seq->at(ref_index) : 'N';
    const bool passes_quality_filter =
        alignment.qualities.modified_qualities.at(query_pos + query_offset) >= _options.min_baseq &&
        alignment.qualities.modified_base_types.at(query_pos + query_offset) >= _options.min_base_type &&
        sequence::IsACGT(query_base);
    // Only update error profile if the base passes quality filters
    const bool is_mismatch = IsMismatch(_options.include_ambiguous_reference_bases, ref_base, query_base);
    if (passes_quality_filter && ref_position_valid) {
      pileup.error_profile.AddAlignedBase(alignment.alignment_metadata);
      if (is_mismatch) {
        pileup.error_profile.AddSubstitution(query_base, ref_base, alignment.alignment_metadata);
      }
    }
    // Update qscore stats by counting both matches and mismatches
    // We don't apply base quality filter for qscore stats because we want to capture all qscores
    const u8 qual = alignment.qualities.modified_qualities.at(query_pos + query_offset);
    pileup.error_profile.AddQscoreData(qual, is_mismatch, alignment.alignment_metadata);
  }
  if (_options.metric_types.has_coverage_metrics) {
    const bool is_concordant_duplex =
        alignment.qualities.modified_base_types.at(query_pos + query_offset) == yc_decode::BaseType::kConcordant &&
        sequence::IsACGT(query_base);
    const bool passes_quality_filter =
        alignment.qualities.modified_qualities.at(query_pos + query_offset) >= _options.min_baseq &&
        alignment.qualities.modified_base_types.at(query_pos + query_offset) >= _options.min_base_type &&
        sequence::IsACGT(query_base);
    pileup.depth_profile.has_coverage = true;
    pileup.depth_profile.AddDepth(is_concordant_duplex, passes_quality_filter);
  }
}

void MetricsWorker::RecordMatchHomopolymerError(
    const Pileup& pileup, Alignment& alignment, const s64 rpos, const u32 query_pos, const u32 query_offset) {
  // Handle homopolymer error if the match is part of a homopolymer and the user has specified to calculate
  // homopolymer metrics
  if (!_options.metric_types.has_accuracy_metrics || !_options.calculate_hp_metrics || !pileup.is_part_of_hp) {
    return;
  }
  const char query_base = io::GetBase(bam_get_seq(alignment.alignment_ptr), query_pos + query_offset);
  const bool ref_position_valid = _ref_seq.has_value() && !_ref_seq->empty() &&
                                  rpos - _current_super_region.start >= 0 &&
                                  rpos - _current_super_region.start < ToSigned(_ref_seq->size());
  const char& ref_base = ref_position_valid ? _ref_seq->at(rpos - _current_super_region.start) : 'N';
  const u8 relative_position_within_hp = pileup.reference_context.relative_position_within_hp;
  alignment.hp_error_profile[rpos - relative_position_within_hp].UpdateMinQuality(
      alignment.qualities.qualities.at(query_pos + query_offset));
  alignment.hp_error_profile[rpos - relative_position_within_hp].UpdateMinBaseType(
      alignment.qualities.base_types.at(query_pos + query_offset));
  const bool passes_hp_quality_filter =
      alignment.qualities.qualities.at(query_pos + query_offset) >= _options.min_baseq &&
      alignment.qualities.base_types.at(query_pos + query_offset) >= _options.min_base_type &&
      sequence::IsACGT(ref_base) && sequence::IsACGT(query_base);
  if (passes_hp_quality_filter && query_base != ref_base) {
    alignment.hp_error_profile[rpos - relative_position_within_hp].AddSubstitution(relative_position_within_hp);
  }
}

RefPosQueryPosPair MetricsWorker::HandleMatchCigar(Alignment& alignment, s64 ref_pos, u32 query_pos, u32 cigar_len) {
  // Process matches and update pileups
  s64 rpos = ref_pos;
  for (u32 j = 0; j < cigar_len; ++j) {
    if (rpos >= _current_super_region.start && rpos < _current_super_region.end) {
      const size_t pileup_head = rpos - _current_super_region.start;
      auto& pileup = _pileups.at(pileup_head);
      if (pileup.valid) {
        RecordMatch(pileup, alignment, rpos, query_pos, j);
        RecordMatchHomopolymerError(pileup, alignment, rpos, query_pos, j);
      }
    }
    ++rpos;
  }
  return {ref_pos + cigar_len, query_pos + cigar_len};
}

void MetricsWorker::RecordInsertion(Pileup& pileup,
                                    const Alignment& alignment,
                                    const std::string& inserted_seq,
                                    const u32 query_pos,
                                    const u32 cigar_len) {
  // We only count an insertion if all bases in the insertion pass the quality filter
  const bool passes_base_quality_filter =
      std::all_of(alignment.qualities.modified_qualities.begin() + query_pos,
                  alignment.qualities.modified_qualities.begin() + query_pos + cigar_len,
                  [this](u8 quality) { return quality >= _options.min_baseq; });
  const bool passes_base_type_filter =
      std::all_of(alignment.qualities.modified_base_types.begin() + query_pos,
                  alignment.qualities.modified_base_types.begin() + query_pos + cigar_len,
                  [this](yc_decode::BaseType base_type) { return base_type >= _options.min_base_type; });
  if (passes_base_quality_filter && passes_base_type_filter) {
    pileup.error_profile.AddInsertion(inserted_seq, alignment.alignment_metadata);
  }
}

void MetricsWorker::RecordInsertionHomopolymerError(const Pileup& pileup,
                                                    Alignment& alignment,
                                                    const s64 ref_pos,
                                                    const std::string& inserted_seq,
                                                    const u32 query_pos,
                                                    const u32 cigar_len) const {
  // Handle homopolymer error if the insertion is part of a homopolymer and the user has specified to calculate
  // homopolymer metrics
  if (!_options.calculate_hp_metrics || !pileup.is_part_of_hp) {
    return;
  }
  const u8 relative_position_within_hp = pileup.reference_context.relative_position_within_hp;
  alignment.hp_error_profile[ref_pos - relative_position_within_hp].AddInsertion(relative_position_within_hp,
                                                                                 inserted_seq);
  const u8 min_unmodified_quality = *std::min_element(alignment.qualities.qualities.begin() + query_pos,
                                                      alignment.qualities.qualities.begin() + query_pos + cigar_len);
  const yc_decode::BaseType min_unmodified_base_type =
      *std::min_element(alignment.qualities.base_types.begin() + query_pos,
                        alignment.qualities.base_types.begin() + query_pos + cigar_len);
  alignment.hp_error_profile[ref_pos - relative_position_within_hp].UpdateMinQuality(min_unmodified_quality);
  alignment.hp_error_profile[ref_pos - relative_position_within_hp].UpdateMinBaseType(min_unmodified_base_type);
}

RefPosQueryPosPair MetricsWorker::HandleInsertionCigar(
    Alignment& alignment, s64 ref_pos, const u32 query_pos, const u32 cigar_len, const bool is_last_cigar_op) {
  // We record an insertion at the nearest mapped base. In most cases, this is the base immediately after
  // the insertion. However, if the insertion is the last CIGAR operation, `rpos` will point to the exclusive
  // end position of the alignment, which is one base after the last mapped base. In this case, we will record
  // the insertion at the last mapped base before the insertion.
  s64 rpos_of_nearest_mapped_base = ref_pos;
  if (is_last_cigar_op && ref_pos > 0) {
    rpos_of_nearest_mapped_base = ref_pos - 1;
  }
  // Check if the position of the nearest mapped base is within the current super region
  if (rpos_of_nearest_mapped_base < _current_super_region.start ||
      rpos_of_nearest_mapped_base >= _current_super_region.end) {
    return {ref_pos, query_pos + cigar_len};
  }
  // Process insertions and update pileups
  if (_options.metric_types.has_accuracy_metrics || _options.metric_types.has_read_metrics) {
    const size_t pileup_head = rpos_of_nearest_mapped_base - _current_super_region.start;
    if (_pileups.at(pileup_head).valid && _options.metric_types.has_accuracy_metrics) {
      auto& pileup = _pileups.at(pileup_head);
      const std::string inserted_seq = io::GetSequence(bam_get_seq(alignment.alignment_ptr), query_pos, cigar_len);
      RecordInsertion(pileup, alignment, inserted_seq, query_pos, cigar_len);
      RecordInsertionHomopolymerError(pileup, alignment, ref_pos, inserted_seq, query_pos, cigar_len);
    }
  }

  return {ref_pos, query_pos + cigar_len};
}

void MetricsWorker::RecordDeletion(
    Pileup& pileup, const Alignment& alignment, const s64 ref_pos, const u32 query_pos, const u32 cigar_len) {
  // We only count a deletion if the bases immediately before and after the deletion pass the quality filter
  const u8 left_base_quality =
      query_pos > 0 ? alignment.qualities.modified_qualities.at(query_pos - 1) : _options.min_baseq;
  const u8 right_base_quality = query_pos < alignment.qualities.modified_qualities.size()
                                    ? alignment.qualities.modified_qualities.at(query_pos)
                                    : _options.min_baseq;
  const yc_decode::BaseType left_base_type =
      query_pos > 0 ? alignment.qualities.modified_base_types.at(query_pos - 1) : _options.min_base_type;
  const yc_decode::BaseType right_base_type = query_pos < alignment.qualities.modified_base_types.size()
                                                  ? alignment.qualities.modified_base_types.at(query_pos)
                                                  : _options.min_base_type;
  if (left_base_quality >= _options.min_baseq && right_base_quality >= _options.min_baseq &&
      left_base_type >= _options.min_base_type && right_base_type >= _options.min_base_type) {
    const bool deletion_within_range = _ref_seq.has_value() && !_ref_seq->empty() &&
                                       ref_pos - _current_super_region.start >= 0 &&
                                       ref_pos - _current_super_region.start + cigar_len < ToSigned(_ref_seq->size());
    const std::string deleted_seq(deletion_within_range && _ref_seq
                                      ? _ref_seq->substr(ref_pos - _current_super_region.start, cigar_len)
                                      : std::string(cigar_len, 'N'));
    pileup.error_profile.AddDeletion(deleted_seq, alignment.alignment_metadata);
  }
}

void MetricsWorker::RecordDeletionHomopolymerError(const Pileup& pileup,
                                                   Alignment& alignment,
                                                   const s64 ref_pos,
                                                   const u32 query_pos) {
  // Handle homopolymer error if the deletion is part of a homopolymer and the user has specified to calculate
  // homopolymer metrics
  if (!_options.calculate_hp_metrics || !pileup.is_part_of_hp) {
    return;
  }
  const u8 relative_position_within_hp = pileup.reference_context.relative_position_within_hp;
  const u8 unmodified_left_base_quality =
      query_pos > 0 ? alignment.qualities.qualities.at(query_pos - 1) : _options.min_baseq;
  const u8 unmodified_right_base_quality = query_pos < alignment.qualities.qualities.size()
                                               ? alignment.qualities.qualities.at(query_pos)
                                               : _options.min_baseq;
  const yc_decode::BaseType& unmodified_left_base_type =
      query_pos > 0 ? alignment.qualities.base_types.at(query_pos - 1) : _options.min_base_type;
  const yc_decode::BaseType& unmodified_right_base_type = query_pos < alignment.qualities.base_types.size()
                                                              ? alignment.qualities.base_types.at(query_pos)
                                                              : _options.min_base_type;
  if (unmodified_left_base_quality >= _options.min_baseq && unmodified_right_base_quality >= _options.min_baseq &&
      unmodified_left_base_type >= _options.min_base_type && unmodified_right_base_type >= _options.min_base_type) {
    alignment.hp_error_profile[ref_pos - relative_position_within_hp].AddDeletion(relative_position_within_hp);
  }
}

RefPosQueryPosPair MetricsWorker::HandleDeletionCigar(Alignment& alignment, s64 ref_pos, u32 query_pos, u32 cigar_len) {
  // Process deletions and update pileups
  if (_options.metric_types.has_accuracy_metrics || _options.metric_types.has_coverage_metrics) {
    if (ref_pos < _current_super_region.start || ref_pos >= _current_super_region.end) {
      return {ref_pos + cigar_len, query_pos};
    }
    const size_t pileup_head = ref_pos - _current_super_region.start;
    if (_pileups.at(pileup_head).valid) {
      auto& pileup = _pileups.at(pileup_head);
      // Deleted positions contribute to coverage metrics even when the user has specified to exclude empty positions
      // because they represent positions that are covered by the read but have no bases aligned to them due to the
      // deletion.
      pileup.depth_profile.has_coverage = true;
      if (_options.metric_types.has_accuracy_metrics) {
        RecordDeletion(pileup, alignment, ref_pos, query_pos, cigar_len);
        RecordDeletionHomopolymerError(pileup, alignment, ref_pos, query_pos);
      }
    }
  }
  return {ref_pos + cigar_len, query_pos};
}

void MetricsWorker::HandleHomopolymerErrors(const Alignment& alignment) {
  for (const auto& [ref_pos, hp_profile] : alignment.hp_error_profile) {
    if (hp_profile.valid && ref_pos >= _current_super_region.start && ref_pos < _current_super_region.end) {
      // Aggregate the homopolymer error profile for this position into the corresponding pileup
      const size_t pileup_head = ref_pos - _current_super_region.start;
      auto& pileup = _pileups.at(pileup_head);
      if (!pileup.is_part_of_hp) {
        continue;
      }
      if (hp_profile.GetMinBaseType() == yc_decode::BaseType::kDiscordant || hp_profile.GetMinQuality() == 0) {
        ++pileup.hp_error_profile.discordant_reads;
      }
      if (hp_profile.GetMinBaseType() < _options.min_base_type || hp_profile.GetMinQuality() < _options.min_baseq ||
          hp_profile.GetMinQuality() == kHomopolymerUninitializedBaseQuality) {
        ++pileup.hp_error_profile.low_quality_reads;
        continue;
      }
      const bool all_insertions_homogeneous = hp_profile.AreInsertionsHomogeneous(pileup.reference_context.hp_base);
      if (!_options.hp_allow_heterogeneous_insertions && !all_insertions_homogeneous) {
        continue;
      }
      ++pileup.hp_error_profile.effective_reads;
      pileup.hp_error_profile.AddErrorCounts(hp_profile);
    }
  }
}

}  // namespace xoos::alignment_metrics
