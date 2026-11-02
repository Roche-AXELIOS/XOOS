#include "metrics/simplex-metrics.h"

#include <fmt/format.h>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/util/hash.h>

#include "core/demux-and-trim-pipeline.h"
#include "metric-filenames.h"     // IWYU pragma: keep
#include "metrics-constants.h"    // IWYU pragma: keep
#include "metrics-constraints.h"  // IWYU pragma: keep
#include "metrics.h"              // IWYU pragma: keep
#include "sequence/sequence-enums.h"

namespace xoos::demux {

using enum SequenceFound;
namespace fs = std::filesystem;

thread_local concurrent::EnumerableThreadLocal<SimplexMetrics> SimplexMetrics::instance{
    std::make_shared<SimplexMetrics>()};

// 2^4 combinations of sid_5p, umi_5p, umi_3p, sid_3p
constexpr u32 kNumSequenceFound = 16;

static vec<SequencesFound> GenerateAllSequencesFound() {
  vec<SequencesFound> result;
  result.reserve(kNumSequenceFound);
  for (const auto sid_5p : {kYes, kNo}) {
    for (const auto sid_3p : {kYes, kNo}) {
      for (const auto umi_5p : {kYes, kNo}) {
        for (const auto umi_3p : {kYes, kNo}) {
          result.emplace_back(SequencesFound{sid_5p, umi_5p, umi_3p, sid_3p});
        }
      }
    }
  }
  return result;
}

const vec<SequencesFound> SimplexMetrics::kAllSequencesFound = GenerateAllSequencesFound();

size_t SequencesFoundHash::operator()(const SequencesFound& sf) const {
  return util::hash::Hash(sf.sid_5p, sf.umi_5p, sf.umi_3p, sf.sid_3p);
}

size_t AdapterMetricHash::operator()(const AdapterMetric& am) const {
  return util::hash::Hash(am.sid, SequencesFoundHash()(am.found));
}

size_t SidCollisionHash::operator()(const SidCollision& sc) const {
  return util::hash::Hash(sc.expected_sid, sc.end, sc.collided_sid);
}

SimplexMetrics& SimplexMetrics::Instance() {
  return *instance.Local();
}

SimplexMetrics SimplexMetrics::SumTotal() {
  auto total = SimplexMetrics{};
  std::function func = [&total](const SimplexMetrics& item) { total.Add(item); };
  instance.ForEach(func);
  return total;
}

static std::string Format(const SequencesFound& sequences_found) {
  const auto* sid_5p = sequences_found.sid_5p == kYes ? "s5" : "xX";
  const auto* sid_3p = sequences_found.sid_3p == kYes ? "s3" : "xX";
  const auto* umi_5p = sequences_found.umi_5p == kYes ? "u5" : "xX";
  const auto* umi_3p = sequences_found.umi_3p == kYes ? "u3" : "xX";
  return fmt::format("{}{}{}{}", sid_5p, sid_3p, umi_5p, umi_3p);
}

SimplexMetrics::SimplexMetrics() {
  for (const auto& sf : kAllSequencesFound) {
    for (size_t sid = 0; sid < kSidPoolSizeHint; ++sid) {
      _adapter_counts[{sid, sf}] = 0;
    }
    _adapter_counts[{std::nullopt, sf}] = 0;
  }

  const auto per_sample_histogram = LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0);
  _sample_untrimmed_simplex_passing_read_length_distr =
      vec(metrics_constraints::max_sid_id_index + 1, per_sample_histogram);
  _sample_untrimmed_full_simplex_passing_read_length_distr =
      vec(metrics_constraints::max_sid_id_index + 1, per_sample_histogram);
  _sample_untrimmed_partial_simplex_passing_read_length_distr =
      vec(metrics_constraints::max_sid_id_index + 1, per_sample_histogram);
  _sample_trimmed_simplex_passing_read_length_distr =
      vec(metrics_constraints::max_sid_id_index + 1, per_sample_histogram);
  _sample_trimmed_full_simplex_passing_read_length_distr =
      vec(metrics_constraints::max_sid_id_index + 1, per_sample_histogram);
  _sample_trimmed_partial_simplex_passing_read_length_distr =
      vec(metrics_constraints::max_sid_id_index + 1, per_sample_histogram);
  _sample_endadapter_position_distr = vec(metrics_constraints::max_sid_id_index + 1, per_sample_histogram);
}

// TODO: replace the std::function with a template argument
u64 SimplexMetrics::TotalCountForPool(const SidPool& sid_pool, const std::function<u64(u32)>& func) {
  auto count_all = u64{0};
  for (const auto& sid : sid_pool) {
    count_all += func(sid.id);
  }
  return count_all;
}

u64 SimplexMetrics::TotalAssignedReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool, [this](const u32 id) { return GetSidReadCount(id); });
}

u64 SimplexMetrics::TotalFoundSids(const SidPool& sid_pool) const {
  // number of found sids with at least one read assigned
  // iterate over _sid_pool and find sids with read_count > 0
  auto num_assigned_sids = u64{0};
  for (const auto& sid : sid_pool) {
    if (GetSidReadCount(sid.id) > 0) {
      num_assigned_sids += 1;
    }
  }
  return num_assigned_sids;
}

u64 SimplexMetrics::TotalFullReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool, [this](const u32 id) { return GetSidFullReadCount(id); });
}

u64 SimplexMetrics::TotalPartialReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool, [this](const u32 id) { return GetSidPartialReadCount(id); });
}

u64 SimplexMetrics::TotalPerfectIndexReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool, [this](const u32 id) { return GetSidPerfectIndexReadCount(id); });
}

u64 SimplexMetrics::TotalSidCount(const SidPool& sid_pool) {
  return sid_pool.size();
}

u64 SimplexMetrics::TotalBothSidDetectedReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool, [this](const u32 id) { return GetSidBothSidDetectedReadCount(id); });
}

u64 SimplexMetrics::TotalSidDiscordantReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool, [this](const u32 id) { return GetSidSidDiscordantReadCount(id); });
}

u64 SimplexMetrics::TotalSidConcordantReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool, [this](const u32 id) { return GetSidSidConcordantReadCount(id); });
}

u64 SimplexMetrics::TotalInputReadCount() const {
  return _total_input_read_count;
}

u64 SimplexMetrics::TotalUnassignedReadCount() const {
  return unassigned_read_count;
}

u64 SimplexMetrics::TotalFilteredReadCount() const {
  return filtered_read_count;
}

u64 SimplexMetrics::GetSidReadCount(const u32 id) const {
  return GetCount(_sid_read_counts, id);
}

u64 SimplexMetrics::GetSidFullReadCount(const u32 id) const {
  return GetCount(_sid_full_read_counts, id);
}

u64 SimplexMetrics::GetSidPartialReadCount(const u32 id) const {
  return GetCount(_sid_partial_read_counts, id);
}

u64 SimplexMetrics::GetSidPassingReadCount(const u32 id) const {
  return GetCount(_sid_passing_read_counts, id);
}

void SimplexMetrics::IncrementSidPassingReadCount(const u32 id) {
  // Drives both the run-level passing count and the per-SID passing count from a single
  // point so their sum can never diverge.
  ++passing_read_count;
  IncrementCount(&_sid_passing_read_counts, id, u64{1});
}

void SimplexMetrics::IncrementSidUmiFailedReadCount(const u32 id) {
  // Drives both the run-level umi_failed count and the per-SID umi_failed count from a single
  // point so their sum can never diverge.
  ++umi_failed_read_count;
  IncrementCount(&_sid_umi_failed_read_counts, id, u64{1});
}

u64 SimplexMetrics::GetSidUmiFailedReadCount(const u32 id) const {
  return GetCount(_sid_umi_failed_read_counts, id);
}

void SimplexMetrics::AddSidPassingBases(const u32 id, const u64 total_bases, const u64 trimmed_bases) {
  // Drives both the run-level and per-SID base totals from a single point so their sums cannot diverge.
  total_base_count += total_bases;
  trimmed_base_count += trimmed_bases;
  IncrementCount(&_sid_total_base_counts, id, total_bases);
  IncrementCount(&_sid_trimmed_base_counts, id, trimmed_bases);
}

u64 SimplexMetrics::GetSidTotalBaseCount(const u32 id) const {
  return GetCount(_sid_total_base_counts, id);
}

u64 SimplexMetrics::GetSidTrimmedBaseCount(const u32 id) const {
  return GetCount(_sid_trimmed_base_counts, id);
}

u64 SimplexMetrics::GetSidPerfectIndexReadCount(const u32 id) const {
  return GetCount(_sid_perfect_index_read_counts, id);
}

u64 SimplexMetrics::GetSidBothSidDetectedReadCount(const u32 id) const {
  return GetCount(_sid_both_sid_detected_read_counts, id);
}

u64 SimplexMetrics::GetSidSidDiscordantReadCount(const u32 id) const {
  return GetCount(_sid_sid_discordant_read_counts, id);
}

u64 SimplexMetrics::GetSidSidConcordantReadCount(const u32 id) const {
  return GetCount(_sid_sid_concordant_read_counts, id);
}

void SimplexMetrics::Add(const SimplexMetrics& other) {
  preassign_passing_read_count += other.preassign_passing_read_count;
  passing_read_count += other.passing_read_count;
  _total_input_read_count += other._total_input_read_count;
  too_short_read_count += other.too_short_read_count;
  too_long_read_count += other.too_long_read_count;
  filtered_read_count += other.filtered_read_count;
  unassigned_read_count += other.unassigned_read_count;
  too_short_trimmed_read_count += other.too_short_trimmed_read_count;
  sid_discordant_discarded_read_count += other.sid_discordant_discarded_read_count;
  full_only_discarded_read_count += other.full_only_discarded_read_count;
  umi_failed_read_count += other.umi_failed_read_count;
  raw_base_count += other.raw_base_count;
  filtered_base_count += other.filtered_base_count;
  unassigned_base_count += other.unassigned_base_count;
  failed_assigned_base_count += other.failed_assigned_base_count;
  total_base_count += other.total_base_count;
  trimmed_base_count += other.trimmed_base_count;

  for (decltype(other._sid_read_counts)::size_type i = 0; i < other._sid_read_counts.size(); ++i) {
    IncrementCount(&_sid_read_counts, i, GetCount(other._sid_read_counts, i));
    IncrementCount(&_sid_full_read_counts, i, GetCount(other._sid_full_read_counts, i));
    IncrementCount(&_sid_partial_read_counts, i, GetCount(other._sid_partial_read_counts, i));
    IncrementCount(&_sid_passing_read_counts, i, GetCount(other._sid_passing_read_counts, i));
    IncrementCount(&_sid_total_base_counts, i, GetCount(other._sid_total_base_counts, i));
    IncrementCount(&_sid_trimmed_base_counts, i, GetCount(other._sid_trimmed_base_counts, i));

    IncrementCount(&_sid_perfect_index_read_counts, i, GetCount(other._sid_perfect_index_read_counts, i));
    IncrementCount(&_sid_both_sid_detected_read_counts, i, GetCount(other._sid_both_sid_detected_read_counts, i));
    IncrementCount(&_sid_sid_discordant_read_counts, i, GetCount(other._sid_sid_discordant_read_counts, i));
    IncrementCount(&_sid_sid_concordant_read_counts, i, GetCount(other._sid_sid_concordant_read_counts, i));
    IncrementCount(&_sid_umi_failed_read_counts, i, GetCount(other._sid_umi_failed_read_counts, i));
  }
  for (const auto& [adapter_metric, _] : other._adapter_counts) {
    IncrementAdapterCount(adapter_metric, other.GetAdapterCount(adapter_metric));
  }

  _total_read_length_distr += other._total_read_length_distr;
  _raw_failed_read_length_distr += other._raw_failed_read_length_distr;
  _untrimmed_passing_read_len_dist += other._untrimmed_passing_read_len_dist;
  _untrimmed_full_passing_read_len_dist += other._untrimmed_full_passing_read_len_dist;
  _untrimmed_partial_passing_read_len_dist += other._untrimmed_partial_passing_read_len_dist;
  _trimmed_passing_read_len_dist += other._trimmed_passing_read_len_dist;
  _trimmed_full_passing_read_len_dist += other._trimmed_full_passing_read_len_dist;
  _trimmed_partial_passing_read_len_dist += other._trimmed_partial_passing_read_len_dist;
  _total_endadapter_position_distr += other._total_endadapter_position_distr;
  for (u32 sid = 0; sid <= metrics_constraints::max_sid_id_index; ++sid) {
    _sample_untrimmed_simplex_passing_read_length_distr[sid] +=
        other._sample_untrimmed_simplex_passing_read_length_distr[sid];
    _sample_untrimmed_full_simplex_passing_read_length_distr[sid] +=
        other._sample_untrimmed_full_simplex_passing_read_length_distr[sid];
    _sample_untrimmed_partial_simplex_passing_read_length_distr[sid] +=
        other._sample_untrimmed_partial_simplex_passing_read_length_distr[sid];
    _sample_trimmed_simplex_passing_read_length_distr[sid] +=
        other._sample_trimmed_simplex_passing_read_length_distr[sid];
    _sample_trimmed_full_simplex_passing_read_length_distr[sid] +=
        other._sample_trimmed_full_simplex_passing_read_length_distr[sid];
    _sample_trimmed_partial_simplex_passing_read_length_distr[sid] +=
        other._sample_trimmed_partial_simplex_passing_read_length_distr[sid];
    _sample_endadapter_position_distr[sid] += other._sample_endadapter_position_distr[sid];
  }

  for (const auto& [sc, count] : other._sid_collisions) {
    _sid_collisions[sc] += count;
  }
}

void SimplexMetrics::IncrementInputReadCount() {
  _total_input_read_count += 1;
}

void SimplexMetrics::RecordInputReadLength(const u64 raw_length) {
  _total_read_length_distr.AddCountToHistogram(raw_length, 1);
}

void SimplexMetrics::RecordRawFailedReadLength(const u64 raw_length) {
  _raw_failed_read_length_distr.AddCountToHistogram(raw_length, 1);
}

void SimplexMetrics::WriteMetrics(const DemuxAndTrimParam& param, const SidPool& sid_pool) const {
  WriteRunMetrics(sid_pool, param);
  WriteSampleMetrics(sid_pool, param);
  WriteSampleAssignmentMetrics(sid_pool, param);
  WriteReadLengthDistributions(sid_pool, param);
}

void SimplexMetrics::WriteReadLengthDistributions(const SidPool& sid_pool, const DemuxAndTrimParam& param) const {
  const auto out_dir = param.out_dir / kMetricsDirectory / kReadLengthDistributionDirectory;
  fs::create_directories(out_dir);
  const io::TsvMetadata tsv_metadata{.command_line_info = param.command_line_info};

  const auto write_single = [&tsv_metadata](
                                const LengthHistogram& histogram, const std::string_view name, const fs::path& path) {
    const histogram::Histograms<u64> histograms{{std::string(name), histogram}};
    histogram::WriteHistogramsToTsv<u64>(
        histograms, path, "read_length", histogram::HistogramBinOutput::kMaxLastBin, {}, tsv_metadata);
  };

  struct PassingHistograms {
    const LengthHistogram& total;
    const LengthHistogram& full;
    const LengthHistogram& partial;
  };

  const auto write_passing = [&tsv_metadata](const PassingHistograms& passing_histograms, const fs::path& path) {
    const histogram::Histograms<u64> histograms{{"passing_reads", passing_histograms.total},
                                                {"full_passing_reads", passing_histograms.full},
                                                {"partial_passing_reads", passing_histograms.partial}};
    histogram::WriteHistogramsToTsv<u64>(
        histograms, path, "read_length", histogram::HistogramBinOutput::kMaxLastBin, {}, tsv_metadata);
  };

  write_single(_total_read_length_distr, "total_reads", out_dir / kTotalReadLengthDistr);
  write_single(_raw_failed_read_length_distr, "raw_failed_reads", out_dir / kRawFailedReadLengthDistr);
  write_passing({.total = _untrimmed_passing_read_len_dist,
                 .full = _untrimmed_full_passing_read_len_dist,
                 .partial = _untrimmed_partial_passing_read_len_dist},
                out_dir / kTotalUntrimmedSimplexPassingReadLengthDistr);
  write_passing({.total = _trimmed_passing_read_len_dist,
                 .full = _trimmed_full_passing_read_len_dist,
                 .partial = _trimmed_partial_passing_read_len_dist},
                out_dir / kTotalTrimmedSimplexPassingReadLengthDistr);
  write_single(_total_endadapter_position_distr, "endadapter_position", out_dir / kTotalEndadapterPositionDistr);

  for (const auto& sid : sid_pool) {
    const auto sample_out_dir = out_dir / SampleReadLengthDistributionDirectory(sid.name);
    fs::create_directories(sample_out_dir);
    write_passing({.total = _sample_untrimmed_simplex_passing_read_length_distr.at(sid.id),
                   .full = _sample_untrimmed_full_simplex_passing_read_length_distr.at(sid.id),
                   .partial = _sample_untrimmed_partial_simplex_passing_read_length_distr.at(sid.id)},
                  sample_out_dir / kSampleUntrimmedSimplexPassingReadLengthDistr);
    write_passing({.total = _sample_trimmed_simplex_passing_read_length_distr.at(sid.id),
                   .full = _sample_trimmed_full_simplex_passing_read_length_distr.at(sid.id),
                   .partial = _sample_trimmed_partial_simplex_passing_read_length_distr.at(sid.id)},
                  sample_out_dir / kSampleTrimmedSimplexPassingReadLengthDistr);
    write_single(_sample_endadapter_position_distr.at(sid.id),
                 "endadapter_position",
                 sample_out_dir / kSampleEndadapterPositionDistr);
  }
}

void SimplexMetrics::WriteSampleMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const {
  const auto sample_output = param.out_dir / kMetricsDirectory / kSampleMetricsFile;
  auto out = std::ofstream{sample_output, std::ofstream::out | std::ofstream::trunc};
  const auto& info = param.command_line_info;

  auto writer = csv::make_tsv_writer(out);

  io::WriteTsvMetadata(out, info);

  vec<std::string_view> run_metrics = {
      kIndexSequence,
      kAssignedReads,
      kAssignedReadsPercentage,
      kFullSimplexReads,
      kFullSimplexReadsPercentage,
      kPartialSimplexReads,
      kPartialSimplexReadsPercentage,
      kBothSidDetectedReads,
      kBothSidDetectedReadsPercentage,
      kSidConcordantReads,
      kSidConcordantReadsPercentage,
      kSidDiscordantReads,
      kSidDiscordantReadsPercentage,
      kPerfectIndexReads,
      kPerfectIndexReadsPercentage,
      kUmiFailedReads,
      kUmiFailedReadsPercentage,
      kPassingReads,
      kPassingReadsPercentage,
      kPassingReadsPercentageAllSamples,
      kTotalBases,
      kTotalBasesPercentageAllSamples,
      kTrimmedBases,
      kTrimmedBasesPercentage,
      kMeanPassingReadLength,
  };

  vec<std::string> header;
  // +2 for row header and unassigned sample
  header.reserve(sid_pool.size() + 2);
  header.emplace_back(kMetric);

  for (const auto& sid : sid_pool) {
    header.emplace_back(sid.name);
  }
  header.emplace_back("Unassigned");
  writer << header;

  vec<vec<std::string>> metric_values(run_metrics.size());
  for (size_t i = 0; i < run_metrics.size(); ++i) {
    // +2 for row header and unassigned sample
    metric_values[i].reserve(sid_pool.size() + 2);
    metric_values[i].emplace_back(run_metrics[i]);
  }
  // Known sample counts. Rows are appended in the same order as run_metrics; the
  // runtime check after building each row guards against the two lists drifting out of sync.
  for (const auto& sid : sid_pool) {
    const auto sid_total_bases = GetSidTotalBaseCount(sid.id);
    const auto sid_trimmed_bases = GetSidTrimmedBaseCount(sid.id);
    const auto sid_assigned_reads = GetSidReadCount(sid.id);
    const auto sid_full_reads = GetSidFullReadCount(sid.id);
    const auto sid_partial_reads = GetSidPartialReadCount(sid.id);
    const auto sid_both_sid_detected_reads = GetSidBothSidDetectedReadCount(sid.id);
    const auto sid_concordant_reads = GetSidSidConcordantReadCount(sid.id);
    const auto sid_discordant_reads = GetSidSidDiscordantReadCount(sid.id);
    const auto sid_perfect_index_reads = GetSidPerfectIndexReadCount(sid.id);
    const auto sid_umi_failed_reads = GetSidUmiFailedReadCount(sid.id);
    const auto sid_passing_reads = GetSidPassingReadCount(sid.id);
    const auto sid_mean_passing_read_length =
        sid_passing_reads == 0 ? 0.0 : static_cast<f64>(sid_total_bases) / static_cast<f64>(sid_passing_reads);

    // Column values for this SID, in the same column order as run_metrics.
    const vec<std::string> sid_values = {
        // index_sequence
        sid.sequence,
        // assigned_reads
        std::to_string(sid_assigned_reads),
        // assigned_reads_percentage
        FormatPercentage(sid_assigned_reads, TotalInputReadCount()),
        // full_simplex_reads
        std::to_string(sid_full_reads),
        // full_simplex_reads_percentage
        FormatPercentage(sid_full_reads, sid_assigned_reads),
        // partial_simplex_reads
        std::to_string(sid_partial_reads),
        // partial_simplex_reads_percentage
        FormatPercentage(sid_partial_reads, sid_assigned_reads),
        // both_sid_detected_reads
        std::to_string(sid_both_sid_detected_reads),
        // both_sid_detected_reads_percentage
        FormatPercentage(sid_both_sid_detected_reads, sid_assigned_reads),
        // sid_concordant_reads
        std::to_string(sid_concordant_reads),
        // sid_concordant_reads_percentage
        FormatPercentage(sid_concordant_reads, sid_both_sid_detected_reads),
        // sid_discordant_reads
        std::to_string(sid_discordant_reads),
        // sid_discordant_reads_percentage
        FormatPercentage(sid_discordant_reads, sid_both_sid_detected_reads),
        // perfect_index_reads
        std::to_string(sid_perfect_index_reads),
        // perfect_index_reads_percentage
        FormatPercentage(sid_perfect_index_reads, sid_both_sid_detected_reads),
        // umi_failed_reads
        std::to_string(sid_umi_failed_reads),
        // umi_failed_reads_percentage
        FormatPercentage(sid_umi_failed_reads, sid_assigned_reads),
        // passing_reads
        std::to_string(sid_passing_reads),
        // passing_reads_percentage
        FormatPercentage(sid_passing_reads, sid_assigned_reads),
        // passing_reads_percentage_all_samples
        FormatPercentage(sid_passing_reads, passing_read_count),
        // total_bases
        std::to_string(sid_total_bases),
        // total_bases_percentage_all_samples
        FormatPercentage(sid_total_bases, total_base_count),
        // trimmed_bases
        std::to_string(sid_trimmed_bases),
        // trimmed_bases_percentage
        FormatPercentage(sid_trimmed_bases, sid_total_bases + sid_trimmed_bases),
        // mean_passing_read_length
        fmt::format("{:.2f}", sid_mean_passing_read_length),
    };
    if (sid_values.size() != run_metrics.size()) {
      throw error::Error(
          "Sample metric row count {} must match run_metrics count {}.", sid_values.size(), run_metrics.size());
    }
    for (size_t row = 0; row < sid_values.size(); ++row) {
      metric_values[row].emplace_back(sid_values[row]);
    }
  }

  // Unknown sid: unassigned reads have no passing bases attributed to a SID.
  // index_sequence is N/A, assigned_reads is the total unassigned count, every
  // count in between is zero, and mean_passing_read_length is 0.00.
  const auto unassigned_reads = TotalUnassignedReadCount();
  vec<std::string> unassigned_values = {
      "NA",
      std::to_string(unassigned_reads),
      FormatPercentage(unassigned_reads, TotalInputReadCount()),
      "0",
      FormatPercentage(0, unassigned_reads),
      "0",
      FormatPercentage(0, unassigned_reads),
      "0",
      FormatPercentage(0, unassigned_reads),
      "0",
      std::string(kNA),
      "0",
      std::string(kNA),
      "0",
      std::string(kNA),
      "0",
      FormatPercentage(0, unassigned_reads),
      "0",
      FormatPercentage(0, unassigned_reads),
      FormatPercentage(0, passing_read_count),
      "0",
      FormatPercentage(0, total_base_count),
      "0",
      std::string(kNA),
      fmt::format("{:.2f}", 0.0),
  };
  for (size_t row = 0; row < unassigned_values.size(); ++row) {
    metric_values[row].emplace_back(unassigned_values[row]);
  }

  // Write run metrics of each sample
  for (size_t i = 0; i < run_metrics.size(); ++i) {
    writer << metric_values[i];
  }

  // Write adapter metrics of each sample
  for (const auto& sf : kAllSequencesFound) {
    vec<std::string> adapter_metric_values;
    // +2 for row header and unassigned sample
    adapter_metric_values.reserve(sid_pool.size() + 2);
    adapter_metric_values.emplace_back(Format(sf));

    // Output adapter match metrics for each sample in sample sheet
    for (const auto& sid : sid_pool) {
      const u64 count = GetAdapterCount({sid.id, sf});
      adapter_metric_values.emplace_back(std::to_string(count));
    }

    // Sum up adapter match counts for each sid not in the sample sheet and unassigned
    u64 unassigned_count = GetAdapterCount({std::nullopt, sf});
    for (u32 id = 0; std::cmp_less(id, _sid_read_counts.size()); ++id) {
      if (id < sid_pool.size()) {
        continue;
      }
      unassigned_count += GetAdapterCount({id, sf});
    }

    adapter_metric_values.emplace_back(std::to_string(unassigned_count));
    writer << adapter_metric_values;
  }
}

void SimplexMetrics::WriteSampleAssignmentMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const {
  const auto output_sample_assignment_metrics = param.out_dir / kMetricsDirectory / kSampleAssignmentMetrics;
  auto out = std::ofstream{output_sample_assignment_metrics, std::ofstream::out | std::ofstream::trunc};
  const auto& info = param.command_line_info;

  auto writer = csv::make_tsv_writer(out);

  io::WriteTsvMetadata(out, info);

  writer << std::make_tuple("expected_sid",
                            "collision_sid",
                            "5'collision",
                            "3'collision",
                            "total_collisions",
                            "both_sid_detected_reads",
                            "percent_sid_collision");
  for (const auto& expected : sid_pool) {
    for (const auto& collided : sid_pool) {
      auto count_5p = GetSidCollisionCount(ReadEnd::k5p, expected.id, collided.id);
      auto count_3p = GetSidCollisionCount(ReadEnd::k3p, expected.id, collided.id);
      auto total_collisions = count_5p + count_3p;
      auto total_reads = GetSidBothSidDetectedReadCount(expected.id);
      auto percentage = FormatPercentage(total_collisions, total_reads);

      writer << std::make_tuple(
          expected.name, collided.name, count_5p, count_3p, total_collisions, total_reads, percentage);
    }
  }
}

void SimplexMetrics::WriteRunMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const {
  const auto run_output = param.out_dir / kMetricsDirectory / kRunMetricsFile;
  auto out = std::ofstream{run_output, std::ofstream::out | std::ofstream::trunc};
  const auto& info = param.command_line_info;

  auto writer = csv::make_tsv_writer(out);

  io::WriteTsvMetadata(out, info);

  auto total_binned_read_count = TotalAssignedReadCount(sid_pool);
  const auto failed_assigned = too_short_trimmed_read_count + sid_discordant_discarded_read_count +
                               full_only_discarded_read_count + umi_failed_read_count;
  // Canonical Key-Value Summary schema: metric_name / value / percentage / denominator.
  // `denominator` names the metric that is the percentage base. Read-count metrics use total_reads;
  // SID-quality metrics use both_sid_detected_reads; base metrics use raw_bases. Reference rows that
  // are trivially 100% of themselves cite their own name. `NA` marks a not-applicable percentage.
  writer << std::make_tuple(kMetricName, kValue, kPercentage, kDenominator);
  writer << std::make_tuple(kTotalReads,
                            _total_input_read_count,
                            FormatPercentage(_total_input_read_count, TotalInputReadCount()),
                            kTotalReads);

  writer << std::make_tuple(kPreassignmentPassingReads,
                            preassign_passing_read_count,
                            FormatPercentage(preassign_passing_read_count, TotalInputReadCount()),
                            kTotalReads);

  writer << std::make_tuple(kFailedReads,
                            TotalFilteredReadCount(),
                            FormatPercentage(TotalFilteredReadCount(), TotalInputReadCount()),
                            kTotalReads);
  writer << std::make_tuple(
      kTooShortReads, too_short_read_count, FormatPercentage(too_short_read_count, TotalInputReadCount()), kTotalReads);
  writer << std::make_tuple(
      kTooLongReads, too_long_read_count, FormatPercentage(too_long_read_count, TotalInputReadCount()), kTotalReads);
  writer << std::make_tuple(kAssignedReads,
                            total_binned_read_count,
                            FormatPercentage(TotalAssignedReadCount(sid_pool), TotalInputReadCount()),
                            kTotalReads);
  writer << std::make_tuple(
      kPassingReads, passing_read_count, FormatPercentage(passing_read_count, TotalInputReadCount()), kTotalReads);
  writer << std::make_tuple(
      kFailedAssignedReads, failed_assigned, FormatPercentage(failed_assigned, TotalInputReadCount()), kTotalReads);
  writer << std::make_tuple(kTooShortTrimmedReads,
                            too_short_trimmed_read_count,
                            FormatPercentage(too_short_trimmed_read_count, TotalInputReadCount()),
                            kTotalReads);
  writer << std::make_tuple(kSidDiscordantDiscardedReads,
                            sid_discordant_discarded_read_count,
                            FormatPercentage(sid_discordant_discarded_read_count, TotalInputReadCount()),
                            kTotalReads);
  writer << std::make_tuple(kFullOnlyDiscardedReads,
                            full_only_discarded_read_count,
                            FormatPercentage(full_only_discarded_read_count, TotalInputReadCount()),
                            kTotalReads);
  writer << std::make_tuple(kUmiFailedReads,
                            umi_failed_read_count,
                            FormatPercentage(umi_failed_read_count, TotalInputReadCount()),
                            kTotalReads);
  writer << std::make_tuple(kUnassignedReads,
                            TotalUnassignedReadCount(),
                            FormatPercentage(TotalUnassignedReadCount(), TotalInputReadCount()),
                            kTotalReads);

  writer << std::make_tuple(kFullSimplexReads,
                            TotalFullReadCount(sid_pool),
                            FormatPercentage(TotalFullReadCount(sid_pool), TotalInputReadCount()),
                            kTotalReads);

  writer << std::make_tuple(kPartialSimplexReads,
                            TotalPartialReadCount(sid_pool),
                            FormatPercentage(TotalPartialReadCount(sid_pool), TotalInputReadCount()),
                            kTotalReads);

  writer << std::make_tuple(kBothSidDetectedReads,
                            TotalBothSidDetectedReadCount(sid_pool),
                            FormatPercentage(TotalBothSidDetectedReadCount(sid_pool), TotalInputReadCount()),
                            kTotalReads);

  writer << std::make_tuple(
      kSidConcordantReads,
      TotalSidConcordantReadCount(sid_pool),
      FormatPercentage(TotalSidConcordantReadCount(sid_pool), TotalBothSidDetectedReadCount(sid_pool)),
      kBothSidDetectedReads);

  writer << std::make_tuple(
      kSidDiscordantReads,
      TotalSidDiscordantReadCount(sid_pool),
      FormatPercentage(TotalSidDiscordantReadCount(sid_pool), TotalBothSidDetectedReadCount(sid_pool)),
      kBothSidDetectedReads);

  writer << std::make_tuple(
      kPerfectIndexReads,
      TotalPerfectIndexReadCount(sid_pool),
      FormatPercentage(TotalPerfectIndexReadCount(sid_pool), TotalBothSidDetectedReadCount(sid_pool)),
      kBothSidDetectedReads);

  // Base counts. raw_bases (all input reads) is the percentage denominator for every base metric.
  // Invariant: raw_bases = total_bases + trimmed_bases + unassigned_bases + failed_assigned_bases + filtered_bases.
  writer << std::make_tuple(kRawBases, raw_base_count, FormatPercentage(raw_base_count, raw_base_count), kRawBases);
  writer << std::make_tuple(
      kTotalBases, total_base_count, FormatPercentage(total_base_count, raw_base_count), kRawBases);
  writer << std::make_tuple(
      kTrimmedBases, trimmed_base_count, FormatPercentage(trimmed_base_count, raw_base_count), kRawBases);
  writer << std::make_tuple(
      kFilteredBases, filtered_base_count, FormatPercentage(filtered_base_count, raw_base_count), kRawBases);
  writer << std::make_tuple(
      kUnassignedBases, unassigned_base_count, FormatPercentage(unassigned_base_count, raw_base_count), kRawBases);
  writer << std::make_tuple(kFailedAssignedBases,
                            failed_assigned_base_count,
                            FormatPercentage(failed_assigned_base_count, raw_base_count),
                            kRawBases);

  // mean_passing_read_length = total_bases / passing_reads (a length, not a percentage), so the
  // percentage and denominator columns are not applicable.
  const auto mean_passing_read_length =
      passing_read_count == 0 ? 0.0 : static_cast<f64>(total_base_count) / static_cast<f64>(passing_read_count);
  writer << std::make_tuple(kMeanPassingReadLength, fmt::format("{:.2f}", mean_passing_read_length), kNA, kNA);

  writer << std::make_tuple(kNumExpectedSids,
                            TotalSidCount(sid_pool),
                            FormatPercentage(TotalSidCount(sid_pool), TotalSidCount(sid_pool)),
                            kNumExpectedSids);
  writer << std::make_tuple(kNumSids,
                            TotalFoundSids(sid_pool),
                            FormatPercentage(TotalFoundSids(sid_pool), TotalSidCount(sid_pool)),
                            kNumExpectedSids);

  for (const auto& sf : kAllSequencesFound) {
    u64 count = 0;
    for (u32 id = 0; id < _sid_read_counts.size(); ++id) {
      count += GetAdapterCount({id, sf});
    }
    count += GetAdapterCount({std::nullopt, sf});
    auto percent = FormatPercentage(count, TotalInputReadCount());
    writer << std::make_tuple(Format(sf), count, percent, kTotalReads);
  }
}

void SimplexMetrics::IncrementAdapterCount(const AdapterMetric& am, const u64 count) {
  _adapter_counts[am] += count;
}

u64 SimplexMetrics::GetAdapterCount(const AdapterMetric& am) const {
  const auto it = _adapter_counts.find(am);
  return it == _adapter_counts.end() ? 0 : it->second;
}

u64 SimplexMetrics::GetSidCollisionCount(const ReadEnd end, const u32 expected_sid, const u32 collided_sid) const {
  const auto it = _sid_collisions.find(SidCollision{end, expected_sid, collided_sid});
  return it == _sid_collisions.end() ? 0 : it->second;
}

}  // namespace xoos::demux
