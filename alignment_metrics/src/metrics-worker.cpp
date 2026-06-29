#include "metrics-worker.h"

#include <algorithm>
#include <memory>
#include <optional>
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

#include "core/alignment.h"
#include "core/pileup.h"
#include "core/super-region.h"
#include "metrics/accuracy-metrics/hp-stat.h"
#include "metrics/metrics.h"
#include "util/homopolymer-util.h"

namespace xoos::alignment_metrics {

BamItrAlignmentProvider::BamItrAlignmentProvider(const SuperRegion& super_region,
                                                 const AlignmentReader& alignment_reader)
    : _alignment_reader(alignment_reader),
      _itr(io::SamItrRegion(
          alignment_reader.idx.get(), alignment_reader.header.get(), CreateHtsRegList(super_region).release(), 1)) {
}

s32 BamItrAlignmentProvider::operator()(bam1_t* const bam1_ptr) const {
  return io::SamItrMultiNext(_alignment_reader.bam.get(), _itr.get(), bam1_ptr);
}

MetricsWorker::MetricsWorker(const AlignmentMetricsOptions& options,
                             const DatasetMetadata& dataset_metadata,
                             const RegionLookupTable& region_lookup_table)
    : _options(options),
      _dataset_metadata(dataset_metadata),
      _region_lookup_table(region_lookup_table),
      _alignment_reader(!_options.bam_input.empty()
                            ? std::make_shared<AlignmentReader>(OpenAlignmentFile(_options.bam_input))
                            : nullptr),
      _metrics(std::make_shared<Metrics>(_options, _dataset_metadata)) {
}

static void InitializeHomopolymerStats(const Homopolymers& hp_regions,
                                       const SuperRegion& super_region,
                                       Metrics& metrics,
                                       std::vector<Pileup>& pileups) {
  for (const auto& [start, end, base] : hp_regions) {
    const auto hp_length = static_cast<u8>(end - start + 1);
    const char hp_base = base;
    // Increment total HP count for this type of HP if calculating HP accuracy metrics
    if (metrics.accuracy_metrics.has_value() && metrics.accuracy_metrics->hp_stats.has_value()) {
      metrics.accuracy_metrics->hp_stats->IncrementHpCount(HpStatKey{hp_base, hp_length});
    }
    for (auto i = start; i <= end; ++i) {
      if (pileups.at(i - ToUnsigned(super_region.start)).valid) {
        pileups.at(i - ToUnsigned(super_region.start)).is_part_of_hp = true;
        if (i == start) {
          pileups.at(i - ToUnsigned(super_region.start)).hp_error_profile.hp_length = hp_length;
          pileups.at(i - ToUnsigned(super_region.start)).hp_error_profile.base = hp_base;
        }
        pileups.at(i - ToUnsigned(super_region.start)).hp_error_profile.relative_position_within_hp =
            static_cast<u8>(i - start);
      }
    }
  }
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
        pileup.valid = false;
        pileup.is_part_of_hp = false;
      }
    }
  }
  _pileups.resize(static_cast<size_t>(super_region.end - super_region.start));
  // Set valid positions in the pileups
  for (const auto& subregion : super_region.subregions) {
    for (s64 i = subregion.start; i < subregion.end; ++i) {
      _pileups.at(i - super_region.start).valid = true;
      // Typically, we skip positions that do not have coverage when calculating coverage metrics
      // but if the user has specified to include empty positions, we will mark them as having coverage in the pileup
      // so that they are included in the coverage histogram calculations.
      if (!_options.exclude_empty_positions) {
        _pileups.at(i - super_region.start).depth_profile.has_coverage = true;
      }
    }
  }
  // Initialize reference sequence for the region if a reference path is provided
  if (ref_seq.has_value()) {
    _ref_seq = std::move(ref_seq);
  }
  // If we are calculating homopolymer metrics, set `is_part_of_hp` to true for all positions in the homopolymer regions
  if (_options.calculate_hp_metrics) {
    const Homopolymers hp_regions = FindReferenceHomopolymers(*_ref_seq,
                                                              _options.min_hp_length,
                                                              _options.max_hp_length,
                                                              super_region,
                                                              _options.hp_subsampling_fraction,
                                                              _options.hp_subsampling_seed);
    InitializeHomopolymerStats(hp_regions, super_region, *_metrics, _pileups);
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
    if (!fs::exists(_options.reference_input)) {
      throw error::Error(fmt::format("Reference file does not exist: {}", _options.reference_input.string()));
    }
    try {
      ref_seq = io::FastaReader(_options.reference_input)
                    .GetSequence(super_region.chromosome,
                                 super_region.start,
                                 super_region.end + static_cast<s64>(_options.reference_padding));
    } catch (const std::exception& e) {
      Logging::Error("Error reading reference sequence for region {}:{}-{}\n{}",
                     super_region.chromosome,
                     super_region.start,
                     super_region.end,
                     e.what());
      throw;
    }
  }
  ProcessRegion(super_region, super_region_id, alignment_provider, std::move(ref_seq));
}

void MetricsWorker::ProcessNoCoverageRegion(const Interval& region) const {
  if (!_options.exclude_empty_positions && _options.metric_types.has_coverage_metrics) {
    _metrics->coverage_metrics->AddEmptyHistogramData(static_cast<u32>(region.end - region.start));
  }
}

void MetricsWorker::ProcessAlignment(Alignment& alignment) {
  // Get the reference position and query position for the first CIGAR operation
  // We extract the relevant trimming information from the CIGAR string
  alignment.Trim(_options.trim_leading_bases, _options.trim_trailing_bases);
  // Check if the alignment should still be processed after trimming
  if (_options.trim_trailing_bases > 0 || _options.trim_leading_bases > 0) {
    // If the alignment no longer spans any reference positions or its length
    // has become zero after trimming, we skip further processing
    if (alignment.reference_length <= 0 || alignment.read_length <= 0) {
      return;
    }
    // Trimming may also cause an alignment to fall outside of regions of interest
    // so we need to check if the trimmed alignment is still in any valid region
    const auto in_valid_region = _region_lookup_table.IsAlignmentInAnyRegion(
        alignment.rpos, alignment.reference_length, _current_super_region_index, _current_subregion_index);
    if (!in_valid_region) {
      return;
    }
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
      alignment.alignment_ptr->core.flag & BAM_FUNMAP;
  // Check if the current alignment is already handled by another subregion in another super region
  // If not we will update read-level metrics for this alignment
  if (_options.metric_types.has_read_metrics && is_safe_to_count) {
    // Process read-level metrics before read-level filtering because we count the total number of reads
    // without filtering for some of the read-level metrics
    _metrics->read_metrics->AddRead(alignment, alignment.alignment_metadata, passes_read_filter);
  }

  if (!passes_read_filter) {
    // Skip further processing if the read does not pass the filter
    return;
  }

  s64 ref_pos = alignment.rpos;
  u32 query_pos = alignment.qpos;

  // Increment total coverage and spanning coverage for homopolymer positions that are covered by this read
  if (_options.metric_types.has_accuracy_metrics && _options.calculate_hp_metrics) {
    const auto ref_end = alignment.rpos + alignment.reference_length;
    for (s64 i = ref_pos; i < alignment.reference_length + alignment.rpos; ++i) {
      const size_t relative_position = i - _current_super_region.start;
      if (relative_position < _pileups.size() && _pileups.at(relative_position).valid &&
          _pileups.at(relative_position).is_part_of_hp &&
          _pileups.at(relative_position).hp_error_profile.relative_position_within_hp == 0) {
        ++_pileups.at(relative_position).hp_error_profile.total_reads;
        // If the read spans the entire HP region with at least `kRequiredAnchorOverlap` bases on each
        // side, then it contributes to the spanning coverage
        const s64 hp_start = i;
        const s64 hp_end = i + _pileups.at(relative_position).hp_error_profile.hp_length;
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

  // Iterate through each CIGAR operation
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
  // After processing all CIGAR operations, handle homopolymer errors for the read if applicable
  if (_options.metric_types.has_accuracy_metrics && _options.calculate_hp_metrics) {
    HandleHomopolymerErrors(alignment);
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
        const char query_base = io::GetBase(bam_get_seq(alignment.alignment_ptr), query_pos + j);
        if (_options.metric_types.has_accuracy_metrics || _options.metric_types.has_coverage_metrics) {
          if (_options.metric_types.has_accuracy_metrics) {
            // Check if rpos is within the region of interest
            const bool ref_position_valid = _ref_seq.has_value() && !_ref_seq->empty() &&
                                            rpos - _current_super_region.start >= 0 &&
                                            rpos - _current_super_region.start < ToSigned(_ref_seq->size());
            const char& ref_base = ref_position_valid ? _ref_seq->at(rpos - _current_super_region.start) : 'N';
            const bool passes_quality_filter =
                alignment.qualities.modified_qualities.at(query_pos + j) >= _options.min_baseq &&
                alignment.qualities.modified_base_types.at(query_pos + j) >= _options.min_base_type &&
                sequence::IsACGT(ref_base) && sequence::IsACGT(query_base);
            // Only update error profile if the base passes quality filters
            if (passes_quality_filter) {
              pileup.error_profile.AddAlignedBase(alignment.alignment_metadata);
              if (query_base != ref_base) {
                pileup.error_profile.AddSubstitution(query_base, ref_base, alignment.alignment_metadata);
              }
            }
            // Update qscore stats by counting both matches and mismatches
            // We don't apply base quality filter for qscore stats because we want to capture all qscores
            const u8 qual = alignment.qualities.modified_qualities.at(query_pos + j);
            pileup.error_profile.AddQscoreData(qual, query_base != ref_base, alignment.alignment_metadata);
            // Handle homopolymer error if the match is part of a homopolymer and the user has specified to calculate
            // homopolymer metrics
            if (_options.calculate_hp_metrics && pileup.is_part_of_hp) {
              const u8 relative_position_within_hp = pileup.hp_error_profile.relative_position_within_hp;
              alignment.hp_error_profile[rpos - relative_position_within_hp].UpdateMinQuality(
                  alignment.qualities.qualities.at(query_pos + j));
              alignment.hp_error_profile[rpos - relative_position_within_hp].UpdateMinBaseType(
                  alignment.qualities.base_types.at(query_pos + j));
              const bool passes_hp_quality_filter =
                  alignment.qualities.qualities.at(query_pos + j) >= _options.min_baseq &&
                  alignment.qualities.base_types.at(query_pos + j) >= _options.min_base_type &&
                  sequence::IsACGT(ref_base) && sequence::IsACGT(query_base);
              if (passes_hp_quality_filter && query_base != ref_base) {
                alignment.hp_error_profile[rpos - relative_position_within_hp].AddSubstitution(
                    relative_position_within_hp);
              }
            }
          }
          if (_options.metric_types.has_coverage_metrics) {
            const bool is_concordant_duplex =
                alignment.qualities.modified_base_types.at(query_pos + j) == yc_decode::BaseType::kConcordant &&
                sequence::IsACGT(query_base);
            const bool passes_quality_filter =
                alignment.qualities.modified_qualities.at(query_pos + j) >= _options.min_baseq &&
                alignment.qualities.modified_base_types.at(query_pos + j) >= _options.min_base_type &&
                sequence::IsACGT(query_base);
            pileup.depth_profile.has_coverage = true;
            pileup.depth_profile.AddDepth(is_concordant_duplex, passes_quality_filter);
          }
        }
      }
    }
    ++rpos;
  }
  return {ref_pos + cigar_len, query_pos + cigar_len};
}

RefPosQueryPosPair MetricsWorker::HandleInsertionCigar(
    Alignment& alignment, s64 ref_pos, const u32 query_pos, const u32 cigar_len, const bool is_last_cigar_op) {
  // We record an insertion at the nearest mapped base. In most cases, this is the base immediately after
  // the insertion. However, if the insertion is the last CIGAR operation, `rpos` will point to the exclusive
  // end position of the alignment, which is one base after the last mapped base. In this case, we will record
  // the insertion at the last mapped base before the insertion.
  const s64 rpos_of_nearest_mapped_base = is_last_cigar_op ? (ref_pos > 0 ? ref_pos - 1 : ref_pos) : ref_pos;
  // Check if the position of the nearest mapped base is within the current super region
  if (rpos_of_nearest_mapped_base < _current_super_region.start ||
      rpos_of_nearest_mapped_base >= _current_super_region.end) {
    return {ref_pos, query_pos + cigar_len};
  }
  // Process insertions and update pileups
  if (_options.metric_types.has_accuracy_metrics || _options.metric_types.has_read_metrics) {
    size_t pileup_head = rpos_of_nearest_mapped_base - _current_super_region.start;
    if (_pileups.at(pileup_head).valid) {
      const std::string inserted_seq = io::GetSequence(bam_get_seq(alignment.alignment_ptr), query_pos, cigar_len);
      if (_options.metric_types.has_accuracy_metrics) {
        auto& pileup = _pileups.at(pileup_head);
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
        // Handle homopolymer error if the insertion is part of a homopolymer and the user has specified to calculate
        // homopolymer metrics
        if (_options.calculate_hp_metrics && pileup.is_part_of_hp) {
          const u8 relative_position_within_hp = pileup.hp_error_profile.relative_position_within_hp;
          alignment.hp_error_profile[ref_pos - relative_position_within_hp].AddInsertion(relative_position_within_hp,
                                                                                         inserted_seq);
          const u8 min_unmodified_quality =
              *std::min_element(alignment.qualities.qualities.begin() + query_pos,
                                alignment.qualities.qualities.begin() + query_pos + cigar_len);
          const yc_decode::BaseType min_unmodified_base_type =
              *std::min_element(alignment.qualities.base_types.begin() + query_pos,
                                alignment.qualities.base_types.begin() + query_pos + cigar_len);
          alignment.hp_error_profile[ref_pos - relative_position_within_hp].UpdateMinQuality(min_unmodified_quality);
          alignment.hp_error_profile[ref_pos - relative_position_within_hp].UpdateMinBaseType(min_unmodified_base_type);
        }
      }
    }
  }

  return {ref_pos, query_pos + cigar_len};
}

RefPosQueryPosPair MetricsWorker::HandleDeletionCigar(Alignment& alignment, s64 ref_pos, u32 query_pos, u32 cigar_len) {
  // Process deletions and update pileups
  if (_options.metric_types.has_accuracy_metrics || _options.metric_types.has_coverage_metrics) {
    if (ref_pos < _current_super_region.start || ref_pos >= _current_super_region.end) {
      return {ref_pos + cigar_len, query_pos};
    }
    size_t pileup_head = ref_pos - _current_super_region.start;
    if (_pileups.at(pileup_head).valid) {
      auto& pileup = _pileups.at(pileup_head);
      // Deleted positions contribute to coverage metrics even when the user has specified to exclude empty positions
      // because they represent positions that are covered by the read but have no bases aligned to them due to the
      // deletion.
      pileup.depth_profile.has_coverage = true;
      if (_options.metric_types.has_accuracy_metrics) {
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
          const bool deletion_within_range =
              _ref_seq.has_value() && !_ref_seq->empty() && ref_pos - _current_super_region.start >= 0 &&
              ref_pos - _current_super_region.start + cigar_len < ToSigned(_ref_seq->size());
          std::string deleted_seq(deletion_within_range && _ref_seq
                                      ? _ref_seq->substr(ref_pos - _current_super_region.start, cigar_len)
                                      : std::string(cigar_len, 'N'));
          pileup.error_profile.AddDeletion(deleted_seq, alignment.alignment_metadata);
        }
        // Handle homopolymer error if the deletion is part of a homopolymer and the user has specified to calculate
        // homopolymer metrics
        if (_options.calculate_hp_metrics && pileup.is_part_of_hp) {
          const u8 relative_position_within_hp = pileup.hp_error_profile.relative_position_within_hp;
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
          if (unmodified_left_base_quality >= _options.min_baseq &&
              unmodified_right_base_quality >= _options.min_baseq &&
              unmodified_left_base_type >= _options.min_base_type &&
              unmodified_right_base_type >= _options.min_base_type) {
            alignment.hp_error_profile[ref_pos - relative_position_within_hp].AddDeletion(relative_position_within_hp);
          }
        }
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
      const bool all_insertions_homogeneous = hp_profile.AreInsertionsHomogeneous(pileup.hp_error_profile.base);
      if (!_options.hp_allow_heterogeneous_insertions && !all_insertions_homogeneous) {
        continue;
      }
      ++pileup.hp_error_profile.effective_reads;
      pileup.hp_error_profile.AddErrorCounts(hp_profile);
    }
  }
}

}  // namespace xoos::alignment_metrics
