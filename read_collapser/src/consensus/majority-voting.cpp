#include "consensus/majority-voting.h"

#include <cstddef>

#include <xoos/util/container-functions.h>

#include "consensus/consensus-matrix.h"
#include "consensus/qscore-calculator.h"
#include "core/read-collapser-options.h"

namespace xoos::read_collapser {

// Finds the majority base and its count from a given count vector.
static std::pair<char, u32> GetMajorityBaseAndCount(const BaseCounts& counts) {
  // Find the iterator to the largest element in the counts vector
  const auto it = std::ranges::max_element(counts);
  // If the largest count is 0, no base covers this position
  if (*it == 0) {
    return {kBaseN, 0};
  }
  // Calculate the index (distance) and cast it to the BaseIndex enum type
  const auto index = static_cast<BaseIndex>(std::distance(counts.begin(), it));
  // Return the base character and its count
  return {ToBase(index), *it};
}

// Check if the base in a singleton cluster of size 2 is duplex amd not N/P/Gap
static bool IsBaseDuplexInSingletonCluster(const ConsensusMatrix& matrix, const size_t column_index) {
  if (matrix.GetReadCount() != 2) {
    // Only applies to clusters of size 2
    return false;
  }
  // Check if both reads are duplex (to ensure it's not two simplex reads)
  if (!matrix.IsDuplex(0) || !matrix.IsDuplex(1)) {
    return false;
  }
  const char r1_base = matrix.GetBase(0, column_index);
  return r1_base != kBaseP && r1_base != kBaseN && r1_base != kBaseGap;
}

// Check if the gap base is tied with the non-gap majority base in either strand
static bool IsGapTiedWithNonGapMajorityBaseInEitherStrand(const ColumnBaseCounts& counts,
                                                          const ColumnConsensusResult::ColumnMetrics& metrics) {
  const auto gap_idx = ToSizeTIndex(kBaseGap);

  // Check if the gap count ties the non-gap majority base count in EITHER strand
  const bool fwd_gap_tied =
      (metrics.fwd_majority_base != kBaseGap && counts.fwd_base_counts[gap_idx] == metrics.fwd_majority_count);
  const bool rev_gap_tied =
      (metrics.rev_majority_base != kBaseGap && counts.rev_base_counts[gap_idx] == metrics.rev_majority_count);

  return fwd_gap_tied || rev_gap_tied;
}

// Check if the gap base is tied with the majority base overall
static bool IsGapTiedWithMajorityBase(const ColumnBaseCounts& counts, const ColumnConsensusResult& result) {
  const size_t majority_idx = ToSizeTIndex(result.majority_base);
  const size_t gap_idx = ToSizeTIndex(kBaseGap);

  // check if the final called base is tied with the gap count
  return counts.overall_base_counts[majority_idx] == counts.overall_base_counts[gap_idx];
}

// Determine the appropriate fallback Q-score based on voting flags and model choice
static u8 DetermineFallbackQscore(const ColumnConsensusResult::ColumnVotingFlags& flags,
                                  const bool enable_legacy_qscore_model) {
  u8 fallback_qscore = enable_legacy_qscore_model ? kLegacyLowConfidenceFallbackQscore : kLowConfidenceFallbackQscore;

  // If the base is simplex OR if all supporting duplex reads were discordant/simplex,
  if (flags.deconvolution_enabled && (flags.is_simplex_base || flags.all_duplex_bases_are_discordant_or_simplex)) {
    // AND the base isn't low confidence for other reasons (strand bias, etc.)
    if (!(flags.has_strand_bias || flags.not_majority || flags.is_singleton || flags.is_hp_insertion_tie)) {
      // Q22 for simplex bases (Q18 if using old qscore model)
      fallback_qscore = enable_legacy_qscore_model ? kLegacySimplexFallbackQscore : kSimplexFallbackQscore;
    }
  }

  if (flags.deconvolution_enabled && flags.parent_daughter_discordant) {
    fallback_qscore = enable_legacy_qscore_model ? kLegacyLowConfidenceFallbackQscore : kLowConfidenceFallbackQscore;
  }

  if (flags.is_singleton) {
    // Q22 for singleton bases (Q0 if using old qscore model)
    fallback_qscore = enable_legacy_qscore_model ? kLegacySingletonFallbackQscore : kSingletonFallbackQscore;
  }

  return fallback_qscore;
}

void ColumnBaseCounts::Reset() {
  std::ranges::fill(fwd_base_counts, 0);
  std::ranges::fill(rev_base_counts, 0);
  std::ranges::fill(overall_base_counts, 0);
  std::ranges::fill(r1_base_counts, 0);
  std::ranges::fill(r2_base_counts, 0);
  read_support = fwd_read_support = rev_read_support = 0;
}

ColumnMajorityVotingWorker::ColumnMajorityVotingWorker(const ReadCollapserOptions& options)
    : _gap_majority_ratio(options.consensus_gap_threshold),
      _majority_ratio(options.consensus_threshold),
      _hd_deconvolution_type(options.duplex_library_type),
      _enable_legacy_qscore_model(options.enable_legacy_qscore_model) {
}

void ColumnMajorityVotingWorker::ComputeColumnConsensusResult(const ConsensusMatrix& matrix,
                                                              const size_t column_index,
                                                              ColumnConsensusResult& result) {
  // Step 1: Accumulate counts for the current column
  AccumulateCounts(matrix, column_index);
  // Early exit if no reads cover this position
  if (_counts.read_support == 0) {
    return;
  }
  // Step 2: Determine the majority base and check initial flags (singleton, majority ratio)
  DetermineMajorityAndInitialFlags(matrix, column_index, result);
  // Step 3: Revert the gap call if low confidence or strand bias is detected
  ApplyGapRules(result);
  // Step 4: Set duplex and homopolymer related flags
  SetDuplexAndHomopolymerFlags(matrix, column_index, result);
  // Step 5: Calculate final quality score based on all flags
  CalculateQualityScore(result);
}

void ColumnMajorityVotingWorker::AccumulateCounts(const ConsensusMatrix& matrix, const size_t column_index) {
  // Ensure the counts are clean before processing the new position
  _counts.Reset();
  for (size_t read_index = 0; read_index < matrix.GetReadCount(); ++read_index) {
    const char base = matrix.GetBase(read_index, column_index);
    // Majority base should not be N or P
    // so we only count A, C, G, T, Gap
    // NOTE: Gap counts towards the total depth and participates in majority voting
    // but N and P do not because they represent positions without read support
    // and majority voting should be taken across only the reads that cover the position
    if (base == kBaseN || base == kBaseP) {
      continue;
    }
    // Count the frequency of each base; count forward and reverse separately
    const auto index = ToSizeTIndex(base);
    ++_counts.overall_base_counts[index];
    ++_counts.read_support;
    if (matrix.GetStrand(read_index) == ReadStrand::kFwd) {
      ++_counts.fwd_base_counts[index];
      ++_counts.fwd_read_support;
    } else {
      ++_counts.rev_base_counts[index];
      ++_counts.rev_read_support;
    }
    // Count R1 and R2 bases for duplex deconvolution if applicable
    if (_hd_deconvolution_type == HDDeconvolutionType::kParentDaughter) {
      using enum DuplexStrand;
      if (matrix.GetDuplexStrand(read_index) == kR1 || matrix.GetDuplexStrand(read_index) == kSimplex) {
        ++_counts.r1_base_counts[index];
      } else if (matrix.GetDuplexStrand(read_index) == kR2) {
        ++_counts.r2_base_counts[index];
      }
    }
  }
}

void ColumnMajorityVotingWorker::DetermineMajorityAndInitialFlags(const ConsensusMatrix& matrix,
                                                                  const size_t column_index,
                                                                  ColumnConsensusResult& result) const {
  auto& metrics = result.metrics;
  auto& flags = result.flags;

  // Find and set majority base and count for each strand
  std::tie(metrics.fwd_majority_base, metrics.fwd_majority_count) = GetMajorityBaseAndCount(_counts.fwd_base_counts);
  std::tie(metrics.rev_majority_base, metrics.rev_majority_count) = GetMajorityBaseAndCount(_counts.rev_base_counts);

  // Find and set overall majority base and count
  const auto [overall_majority_base, overall_majority_count] = GetMajorityBaseAndCount(_counts.overall_base_counts);
  result.majority_base = overall_majority_base;

  // Set the total read support
  result.read_support = _counts.read_support;

  // Set metrics needed for optional FASTQ tags
  metrics.fwd_read_support = _counts.fwd_read_support;
  metrics.rev_read_support = _counts.rev_read_support;

  // Check if the position is a singleton (cluster size 1 is not reliable)
  flags.is_singleton = (result.read_support == 1);
  // Check if the majority base meets the majority ratio threshold
  flags.not_majority = (static_cast<f64>(overall_majority_count) / result.read_support <= _majority_ratio);
  // Check for strand bias: both strands have different majorities and neither is N
  flags.has_strand_bias = metrics.fwd_majority_base != kBaseN && metrics.rev_majority_base != kBaseN &&
                          metrics.fwd_majority_base != metrics.rev_majority_base;

  flags.deconvolution_enabled = (_hd_deconvolution_type != HDDeconvolutionType::kNone);
  const bool is_part_of_hp = matrix.IsPartOfHomopolymer(column_index);
  const bool is_adjacent_to_hp =
      (column_index > 0 && matrix.IsPartOfHomopolymer(column_index - 1)) ||
      (column_index < matrix.GetConsensusLength() - 1 && matrix.IsPartOfHomopolymer(column_index + 1));
  flags.apply_hp_quality_override = flags.deconvolution_enabled && (is_part_of_hp || is_adjacent_to_hp);
  // If deconvolution is enabled and the cluster is a singleton cluster of a duplex read (i.e. exactly two reads after
  // deconvolution), then instead of breaking ties by lexicographic order, we always pick R1 as the majority base if
  // R1 and R2 do not agree. In the case of a singleton cluster of a duplex read, this is essentially the same as
  // intramolecular consensus and picking R1 as the consensus base is consistent with the behavior of demux
  if (flags.deconvolution_enabled && IsBaseDuplexInSingletonCluster(matrix, column_index)) {
    const auto r1_base = matrix.GetBase(0, column_index);
    result.majority_base = r1_base;
    // Update metrics to reflect the new majority base
    const auto r1_index = ToSizeTIndex(r1_base);
    metrics.fwd_majority_count = _counts.fwd_base_counts[r1_index];
    metrics.rev_majority_count = _counts.rev_base_counts[r1_index];
    metrics.fwd_majority_base = r1_base;
    metrics.rev_majority_base = r1_base;
  }
}

void ColumnMajorityVotingWorker::ApplyGapRules(ColumnConsensusResult& result) {
  // Return early if the majority base is NOT a gap
  if (result.majority_base != kBaseGap) {
    return;
  }

  auto& metrics = result.metrics;
  auto& flags = result.flags;

  // For parent-parent duplex, we skip the strand bias check if
  // (1) all reads in the cluster are duplex (as indicated by the same number of forward and reverse reads after
  // deconvolution), and (2) the position is either a part of or adjacent to a homopolymer, and (3) the gap base is
  // the majority in one strand and tied with another base in the other strand
  //
  // This is to address cases like
  //
  // Read 1: ATCTATTTT
  // Qual 1: ~~~~~~~~~
  // Read 2: ATCTATTTTT
  // Qual 2: ~~~~~!~~~~
  // After deconvolution, we have:
  // Read 1 R1: ATCTA-TTTT (fwd)
  // Read 1 R2: ATCTA-TTTT (rev)
  // Read 2 R1: ATCTATTTTT (fwd)
  // Read 2 R2: ATCTA-TTTT (rev)
  //
  // The strand bias check would normally prevent us from calling a gap here, but since homopolymers have
  // a tendency to contain insertions but only in one strand of the duplex read, we allow the gap to be called
  if (_hd_deconvolution_type == HDDeconvolutionType::kParentParent &&
      metrics.fwd_read_support == metrics.rev_read_support && flags.apply_hp_quality_override &&
      IsGapTiedWithNonGapMajorityBaseInEitherStrand(_counts, metrics)) {
    flags.has_strand_bias = false;
  }

  // Check the consensus gap threshold
  const auto gap_index = ToSizeTIndex(kBaseGap);
  const auto& base_counts = _counts.overall_base_counts[gap_index];

  flags.not_gap_majority = static_cast<f64>(base_counts) / result.read_support <= _gap_majority_ratio;

  // Final decision: revert gap call if low confidence
  // We only call a gap outside of homopolymers if it is a gap majority AND both strands agree
  // Otherwise, we use the second most frequent base
  if (flags.has_strand_bias || flags.not_gap_majority) {
    const auto second_largest = util::container::FindSecondLargestUnique(_counts.overall_base_counts.begin(),
                                                                         _counts.overall_base_counts.end());
    if (second_largest != _counts.overall_base_counts.end() && *second_largest > 0) {
      // Set the majority base to the second most frequent base
      const auto index = static_cast<size_t>(std::distance(_counts.overall_base_counts.begin(), second_largest));
      const auto new_majority_base = ToBase(static_cast<BaseIndex>(index));
      result.majority_base = new_majority_base;
      flags.not_majority = (static_cast<f64>(*second_largest) / result.read_support <= _majority_ratio);

      // Update metrics to reflect the new majority base
      metrics.fwd_majority_count = _counts.fwd_base_counts[index];
      metrics.rev_majority_count = _counts.rev_base_counts[index];
      metrics.fwd_majority_base = new_majority_base;
      metrics.rev_majority_base = new_majority_base;
    }
  }
}

void ColumnMajorityVotingWorker::SetDuplexAndHomopolymerFlags(const ConsensusMatrix& matrix,
                                                              const size_t column_index,
                                                              ColumnConsensusResult& result) const {
  auto& flags = result.flags;

  // We can reduce HP insertion calls by checking if the majority base is tied with the gap base
  // or if the second most frequent base at the position is a gap when the position is a part of or adjacent to a
  // homopolymer. But in both cases, we do not have enough evidence to confidently call a gap, so we call a base with
  // Q0 to indicate a low confidence insertion in a homopolymer region.
  if (result.majority_base != kBaseGap && flags.apply_hp_quality_override &&
      IsGapTiedWithMajorityBase(_counts, result)) {
    // If the majority base is tied with the gap base, we assume it is a homopolymer insertion
    // and assign Q0. This is because homopolymer insertions are often caused by sequencing errors,
    // but we also cannot confidently call a gap without enough supporting reads. Downstream
    // applications can choose to mask the homopolymer containing a Q0 base.
    flags.is_hp_insertion_tie = true;
  }
  // Check for duplex and simplex status.
  if (flags.deconvolution_enabled) {
    for (size_t j = 0; j + 1 < matrix.GetReadCount(); ++j) {
      if (matrix.GetBase(j, column_index) != kBaseP && matrix.GetBase(j + 1, column_index) != kBaseP &&
          matrix.GetDuplexStrand(j) == DuplexStrand::kR1 && matrix.GetDuplexStrand(j + 1) == DuplexStrand::kR2) {
        // If we find an R1/R2 pair with real base support, it's not a simplex base
        flags.is_simplex_base = false;
        break;
      }
    }
    for (size_t j = 0; j < matrix.GetReadCount(); ++j) {
      if (matrix.IsDuplex(j) && matrix.GetBaseType(j, column_index) == yc_decode::BaseType::kConcordant) {
        // If any base is concordant, then the overall quality override shouldn't apply
        flags.all_duplex_bases_are_discordant_or_simplex = false;
        break;
      }
    }
    if (_hd_deconvolution_type == HDDeconvolutionType::kParentDaughter) {
      // Find R1 and R2 majorities
      const auto [r1_majority, r1_count] = GetMajorityBaseAndCount(_counts.r1_base_counts);
      const auto [r2_majority, r2_count] = GetMajorityBaseAndCount(_counts.r2_base_counts);

      if (r1_majority != kBaseN && r2_majority != kBaseN && r1_majority != r2_majority) {
        // If the majorities of R1 and R2 differ, we assign Q5 (Q0 if using old qscore model) to the base
        // This is to prevent calling a gap in a parent-daughter cluster where R1 and R2 do not agree
        flags.parent_daughter_discordant = true;
      }
    }
  }
}

void ColumnMajorityVotingWorker::CalculateQualityScore(ColumnConsensusResult& result) {
  const auto& metrics = result.metrics;
  const auto& flags = result.flags;

  // Determine the baseline statistical q score
  const auto majority_index = ToSizeTIndex(result.majority_base);
  const auto statistical_qscore =
      _enable_legacy_qscore_model
          ? QscoreCalculator::ComputeQualLegacy(result.read_support, _counts.overall_base_counts[majority_index])
          : ComputeQual(metrics.fwd_read_support,
                        metrics.fwd_majority_count,
                        metrics.rev_read_support,
                        metrics.rev_majority_count);

  // Define the conditions that trigger the low-confidence fallback score
  const bool use_singleton_fallback = flags.is_singleton;
  const bool use_low_confidence_fallback = flags.not_majority || flags.has_strand_bias || flags.not_gap_majority;
  const bool duplex_discordant_condition = flags.is_hp_insertion_tie || flags.is_simplex_base ||
                                           flags.parent_daughter_discordant ||
                                           flags.all_duplex_bases_are_discordant_or_simplex;
  const bool use_duplex_discordant_fallback = flags.deconvolution_enabled && duplex_discordant_condition;

  // Final Q score decision
  if (use_singleton_fallback || use_low_confidence_fallback || use_duplex_discordant_fallback) {
    result.qscore = DetermineFallbackQscore(flags, _enable_legacy_qscore_model);
  } else {
    result.qscore = statistical_qscore;
  }
}

}  // namespace xoos::read_collapser
