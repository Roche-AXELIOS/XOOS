#include "metrics/simplex-metrics.h"

#include <fmt/format.h>

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <memory>

#include <csv.hpp>

#include "core/demux-and-trim-pipeline.h"
#include "metric-filenames.h"
#include "metrics-constraints.h"
#include "metrics.h"

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

SimplexMetrics& SimplexMetrics::Instance() { return *instance.Local(); }

SimplexMetrics SimplexMetrics::SumTotal() {
  auto total = SimplexMetrics{};
  std::function func = [&total](const SimplexMetrics& item) { total.Add(item); };
  instance.ForEach(func);
  return total;
}

static std::string Format(const SequencesFound& sequences_found) {
  auto sid_5p = sequences_found.sid_5p == kYes ? "s5" : "xX";
  auto sid_3p = sequences_found.sid_3p == kYes ? "s3" : "xX";
  auto umi_5p = sequences_found.umi_5p == kYes ? "u5" : "xX";
  auto umi_3p = sequences_found.umi_3p == kYes ? "u3" : "xX";
  return fmt::format("{}{}{}{}", sid_5p, sid_3p, umi_5p, umi_3p);
}

SimplexMetrics::SimplexMetrics() {
  for (const auto& sf : kAllSequencesFound) {
    for (size_t sid = 0; sid < kSidPoolSizeHint; ++sid) {
      _adapter_counts[{sid, sf}] = 0;
    }
    _adapter_counts[{std::nullopt, sf}] = 0;
  }

  _untrimmed_full_read_len_dist = vec(metrics_constraints::max_sid_id_index + 1,
                                      histogram::Histogram<u64>(metrics_constraints::max_logged_read_length + 1, 0));
  _untrimmed_partial_read_len_dist = vec(metrics_constraints::max_sid_id_index + 1,
                                         histogram::Histogram<u64>(metrics_constraints::max_logged_read_length + 1, 0));
  _trimmed_full_read_len_dist = vec(metrics_constraints::max_sid_id_index + 1,
                                    histogram::Histogram<u64>(metrics_constraints::max_logged_read_length + 1, 0));
  _trimmed_partial_read_len_dist = vec(metrics_constraints::max_sid_id_index + 1,
                                       histogram::Histogram<u64>(metrics_constraints::max_logged_read_length + 1, 0));
}

// TODO: replace the std::function with a template argument
u64 SimplexMetrics::TotalCountForPool(const SidPool& sid_pool, const std::function<u64(u32)>& func) {
  auto count_all = u64{0};
  for (const auto& [_, sid] : sid_pool) {
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
  for (const auto& [_, sid] : sid_pool) {
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

u64 SimplexMetrics::TotalIndexHoppingReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool,
                           [this, sid_pool](const u32 id) { return GetSidIndexHoppingReadCount(sid_pool, id); });
}

u64 SimplexMetrics::TotalPerfectIndexReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool, [this](const u32 id) { return GetSidPerfectIndexReadCount(id); });
}

u64 SimplexMetrics::TotalSidCount(const SidPool& sid_pool) { return sid_pool.size(); }

u64 SimplexMetrics::TotalBothSidDetectedReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool, [this](const u32 id) { return GetSidBothSidDetectedReadCount(id); });
}

u64 SimplexMetrics::TotalSidDiscordantReadCount(const SidPool& sid_pool) const {
  return TotalCountForPool(sid_pool, [this](const u32 id) { return GetSidSidDiscordantReadCount(id); });
}

u64 SimplexMetrics::TotalReadCount() const { return _total_read_count; }

u64 SimplexMetrics::TotalInputReadCount() const { return _total_input_read_count; }

u64 SimplexMetrics::TotalUnassignedReadCount(const SidPool& sid_pool) const {
  return TotalReadCount() - TotalAssignedReadCount(sid_pool);
}

u64 SimplexMetrics::TotalFilteredReadCount() const { return TotalInputReadCount() - TotalReadCount(); }

f32 SimplexMetrics::PercentFilteredReads() const {
  return CalculatePercentage(TotalFilteredReadCount(), TotalInputReadCount());
}

f32 SimplexMetrics::PercentAssignedReads(const SidPool& sid_pool) const {
  return CalculatePercentage(TotalAssignedReadCount(sid_pool), TotalInputReadCount());
}

f32 SimplexMetrics::PercentIndexHoppingReads(const SidPool& sid_pool) const {
  return CalculatePercentage(TotalIndexHoppingReadCount(sid_pool), TotalBothSidDetectedReadCount(sid_pool));
}

f32 SimplexMetrics::PercentFullReads(const SidPool& sid_pool) const {
  return CalculatePercentage(TotalFullReadCount(sid_pool), TotalInputReadCount());
}

f32 SimplexMetrics::PercentPartialReads(const SidPool& sid_pool) const {
  return CalculatePercentage(TotalPartialReadCount(sid_pool), TotalInputReadCount());
}

f32 SimplexMetrics::PercentPerfectIndexReads(const SidPool& sid_pool) const {
  return CalculatePercentage(TotalPerfectIndexReadCount(sid_pool), TotalBothSidDetectedReadCount(sid_pool));
}

f32 SimplexMetrics::PercentBothSidDetectedReads(const SidPool& sid_pool) const {
  return CalculatePercentage(TotalBothSidDetectedReadCount(sid_pool), TotalInputReadCount());
}

f32 SimplexMetrics::PercentSidDiscordantReads(const SidPool& sid_pool) const {
  return CalculatePercentage(TotalSidDiscordantReadCount(sid_pool), TotalBothSidDetectedReadCount(sid_pool));
}

u64 SimplexMetrics::GetSidReadCount(const u32 id) const { return GetCount(_sid_read_counts, id); }

u64 SimplexMetrics::GetSidFullReadCount(const u32 id) const { return GetCount(_sid_full_read_counts, id); }

u64 SimplexMetrics::GetSidIndexHoppingReadCount(const SidPool& sid_pool, const u32 id) const {
  u64 count = 0;
  for (const auto& [sid, _] : sid_pool) {
    count += GetSidCollisionCount(ReadEnd::k5p, id, sid);
    count += GetSidCollisionCount(ReadEnd::k3p, id, sid);
  }
  return count;
}

u64 SimplexMetrics::GetSidPartialReadCount(const u32 id) const { return GetCount(_sid_partial_read_counts, id); }

u64 SimplexMetrics::GetSidPerfectIndexReadCount(const u32 id) const {
  return GetCount(_sid_perfect_index_read_counts, id);
}

u64 SimplexMetrics::GetSidBothSidDetectedReadCount(const u32 id) const {
  return GetCount(_sid_both_sid_detected_read_counts, id);
}

u64 SimplexMetrics::GetSidSidDiscordantReadCount(const u32 id) const {
  return GetCount(_sid_sid_discordant_read_counts, id);
}

void SimplexMetrics::Add(const SimplexMetrics& other) {
  _total_read_count += other._total_read_count;
  _total_input_read_count += other._total_input_read_count;
  _min_trimmed_read_len_count += other._min_trimmed_read_len_count;

  for (decltype(other._sid_read_counts)::size_type i = 0; i < other._sid_read_counts.size(); ++i) {
    IncrementCount(&_sid_read_counts, i, GetCount(other._sid_read_counts, i));
    IncrementCount(&_sid_full_read_counts, i, GetCount(other._sid_full_read_counts, i));
    IncrementCount(&_sid_partial_read_counts, i, GetCount(other._sid_partial_read_counts, i));

    IncrementCount(&_sid_perfect_index_read_counts, i, GetCount(other._sid_perfect_index_read_counts, i));
    IncrementCount(&_sid_both_sid_detected_read_counts, i, GetCount(other._sid_both_sid_detected_read_counts, i));
    IncrementCount(&_sid_sid_discordant_read_counts, i, GetCount(other._sid_sid_discordant_read_counts, i));
  }
  for (const auto& [adapter_metric, _] : other._adapter_counts) {
    IncrementAdapterCount(adapter_metric, other.GetAdapterCount(adapter_metric));
  }

  for (u32 sid = 0; sid <= metrics_constraints::max_sid_id_index; ++sid) {
    _untrimmed_full_read_len_dist[sid] += other._untrimmed_full_read_len_dist[sid];
    _untrimmed_partial_read_len_dist[sid] += other._untrimmed_partial_read_len_dist[sid];
    _trimmed_full_read_len_dist[sid] += other._trimmed_full_read_len_dist[sid];
    _trimmed_partial_read_len_dist[sid] += other._trimmed_partial_read_len_dist[sid];
  }

  for (const auto& [sc, count] : other._sid_collisions) {
    _sid_collisions[sc] += count;
  }
}

void SimplexMetrics::IncrementInputReadCount() { _total_input_read_count += 1; }

void SimplexMetrics::IncrementMinTrimmedReadLenFilteredCount() { _min_trimmed_read_len_count += 1; }

void SimplexMetrics::WriteMetrics(const DemuxAndTrimParam& param, const SidPool& sid_pool) const {
  WriteRunMetrics(sid_pool, param);
  WriteSampleMetrics(sid_pool, param);
  WriteSampleAssignmentMetrics(sid_pool, param);
  WriteReadLengthDistributions(sid_pool, param);
}

void SimplexMetrics::WriteReadLengthDistributions(const SidPool& sid_pool, const DemuxAndTrimParam& param) const {
  const auto out_dir = param.out_dir / kMetricsDirectory;
  const auto& comments = param.comments;

  auto build_histograms = [&sid_pool](const vec<histogram::Histogram<u64>>& per_sid_hist) {
    histogram::Histograms<u64> histograms;
    histograms.reserve(sid_pool.size());
    for (const auto& [_, barcode] : sid_pool) {
      histograms.emplace_back(barcode.name, per_sid_hist.at(barcode.id));
    }
    return histograms;
  };

  auto write_one = [&comments, &build_histograms](const vec<histogram::Histogram<u64>>& per_sid_hist,
                                                  const fs::path& path) {
    const auto histograms = build_histograms(per_sid_hist);
    histogram::WriteHistograms<u64>(histograms, path, "read_length", histogram::kMaxLastBin, {}, comments);
  };

  write_one(_untrimmed_full_read_len_dist, out_dir / "untrimmed_full_read_len_dist.tsv");
  write_one(_untrimmed_partial_read_len_dist, out_dir / "untrimmed_partial_read_len_dist.tsv");
  write_one(_trimmed_full_read_len_dist, out_dir / "trimmed_full_read_len_dist.tsv");
  write_one(_trimmed_partial_read_len_dist, out_dir / "trimmed_partial_read_len_dist.tsv");
}

void SimplexMetrics::WriteSampleMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const {
  const auto sample_output = param.out_dir / kMetricsDirectory / kSampleMetricsFile;
  auto out = std::ofstream{sample_output, std::ofstream::out | std::ofstream::trunc};
  const auto& comments = param.comments;

  auto writer = csv::make_tsv_writer(out);

  io::WriteTsvComments(writer, comments);

  vec<std::string> run_metrics = {
      kIndexSequence,        kAssignedReads,      kFullReads,         kPartialReads,
      kBothSidDetectedReads, kSidDiscordantReads, kIndexHoppingReads, kPerfectIndexReads,
  };

  vec<std::string> header;
  // +2 for row header and unassigned sample
  header.reserve(sid_pool.size() + 2);
  header.emplace_back(kMetric);

  for (const auto& [_, sid] : sid_pool) {
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
  // Known sample counts
  for (const auto& [_, sid] : sid_pool) {
    metric_values[0].emplace_back(sid.sequence);
    metric_values[1].emplace_back(std::to_string(GetSidReadCount(sid.id)));
    metric_values[2].emplace_back(std::to_string(GetSidFullReadCount(sid.id)));
    metric_values[3].emplace_back(std::to_string(GetSidPartialReadCount(sid.id)));
    metric_values[4].emplace_back(std::to_string(GetSidBothSidDetectedReadCount(sid.id)));
    metric_values[5].emplace_back(std::to_string(GetSidSidDiscordantReadCount(sid.id)));
    metric_values[6].emplace_back(std::to_string(GetSidIndexHoppingReadCount(sid_pool, sid.id)));
    metric_values[7].emplace_back(std::to_string(GetSidPerfectIndexReadCount(sid.id)));
  }

  // Unknown sid
  metric_values[0].emplace_back("NA");
  metric_values[1].emplace_back(std::to_string(TotalUnassignedReadCount(sid_pool)));
  metric_values[2].emplace_back("0");
  metric_values[3].emplace_back("0");
  metric_values[4].emplace_back("0");
  metric_values[5].emplace_back("0");
  metric_values[6].emplace_back("0");
  metric_values[7].emplace_back("0");

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
    for (const auto& [_, sid] : sid_pool) {
      const u64 count = GetAdapterCount({sid.id, sf});
      adapter_metric_values.emplace_back(std::to_string(count));
    }

    // Sum up adapter match counts for each sid not in the sample sheet and unassigned
    u64 unassigned_count = GetAdapterCount({std::nullopt, sf});
    for (u32 id = 0; std::cmp_less(id, _sid_read_counts.size()); ++id) {
      if (sid_pool.contains(id)) {
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
  const auto& comments = param.comments;

  vec<Barcode> barcodes;
  barcodes.reserve(sid_pool.size());
  std::ranges::transform(sid_pool, std::back_insert_iterator(barcodes), [](const auto& p) { return p.second; });
  std::ranges::stable_sort(barcodes, [](const auto& a, const auto& b) { return a.id < b.id; });

  auto writer = csv::make_tsv_writer(out);

  io::WriteTsvComments(writer, comments);

  writer << std::make_tuple("expected_sid", "collision_sid", "5'collision", "3'collision", "total_collisions",
                            "both_sid_detected_reads", "percent_sid_collision");
  for (const auto& expected : barcodes) {
    for (const auto& collided : barcodes) {
      auto count_5p = GetSidCollisionCount(ReadEnd::k5p, expected.id, collided.id);
      auto count_3p = GetSidCollisionCount(ReadEnd::k3p, expected.id, collided.id);
      auto total_collisions = count_5p + count_3p;
      auto total_reads = GetSidBothSidDetectedReadCount(expected.id);
      auto percentage = FormatPercent(CalculatePercentage(total_collisions, total_reads), 4);

      writer << std::make_tuple(expected.name, collided.name, count_5p, count_3p, total_collisions, total_reads,
                                percentage);
    }
  }
}

std::string SimplexMetrics::FormatPercent(const f32 data, const s32 precision) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << data;
  return ss.str();
}

void SimplexMetrics::WriteRunMetrics(const SidPool& sid_pool, const DemuxAndTrimParam& param) const {
  const auto run_output = param.out_dir / kMetricsDirectory / kRunMetricsFile;
  auto out = std::ofstream{run_output, std::ofstream::out | std::ofstream::trunc};
  const auto& comments = param.comments;

  auto writer = csv::make_tsv_writer(out);

  io::WriteTsvComments(writer, comments);

  auto total_binned_read_count = TotalAssignedReadCount(sid_pool);
  writer << std::make_tuple(kMetricsName, kCount, kPercentage);
  writer << std::make_tuple(kTotalReads, _total_input_read_count,
                            FormatPercent(CalculatePercentage(_total_input_read_count, TotalInputReadCount())));

  writer << std::make_tuple(kPreassignmentPassingReads, _total_read_count,
                            FormatPercent(CalculatePercentage(_total_read_count, TotalInputReadCount())));

  writer << std::make_tuple(kFailedReads, TotalFilteredReadCount(), FormatPercent(PercentFilteredReads()));
  writer << std::make_tuple(kTooShortReads, (TotalFilteredReadCount() - _min_trimmed_read_len_count),
                            FormatPercent(CalculatePercentage(TotalFilteredReadCount() - _min_trimmed_read_len_count,
                                                              TotalInputReadCount())));
  writer << std::make_tuple(kTooShortTrimmedReads, _min_trimmed_read_len_count,
                            FormatPercent(CalculatePercentage(_min_trimmed_read_len_count, TotalInputReadCount())));

  writer << std::make_tuple(kAssignedReads, total_binned_read_count, FormatPercent(PercentAssignedReads(sid_pool)));

  writer << std::make_tuple(
      kUnassignedReads, TotalUnassignedReadCount(sid_pool),
      FormatPercent(CalculatePercentage(TotalUnassignedReadCount(sid_pool), TotalInputReadCount())));

  writer << std::make_tuple(kFullReads, TotalFullReadCount(sid_pool), FormatPercent(PercentFullReads(sid_pool)));

  writer << std::make_tuple(kPartialReads, TotalPartialReadCount(sid_pool),
                            FormatPercent(PercentPartialReads(sid_pool)));

  writer << std::make_tuple(kBothSidDetectedReads, TotalBothSidDetectedReadCount(sid_pool),
                            FormatPercent(PercentBothSidDetectedReads(sid_pool)));

  writer << std::make_tuple(kSidDiscordantReads, TotalSidDiscordantReadCount(sid_pool),
                            FormatPercent(PercentSidDiscordantReads(sid_pool)));

  writer << std::make_tuple(kIndexHoppingReads, TotalIndexHoppingReadCount(sid_pool),
                            FormatPercent(PercentIndexHoppingReads(sid_pool)));

  writer << std::make_tuple(kPerfectIndexReads, TotalPerfectIndexReadCount(sid_pool),
                            FormatPercent(PercentPerfectIndexReads(sid_pool)));

  writer << std::make_tuple(kNumExpectedSids, TotalSidCount(sid_pool), FormatPercent(100.0));
  writer << std::make_tuple(
      kNumSids, TotalFoundSids(sid_pool),
      FormatPercent(100.0f * static_cast<f32>(TotalFoundSids(sid_pool)) / static_cast<f32>(TotalSidCount(sid_pool))));

  for (const auto& sf : kAllSequencesFound) {
    u64 count = 0;
    for (u32 id = 0; id < _sid_read_counts.size(); ++id) {
      count += GetAdapterCount({id, sf});
    }
    count += GetAdapterCount({std::nullopt, sf});
    auto percent = FormatPercent(CalculatePercentage(count, TotalInputReadCount()));
    writer << std::make_tuple(Format(sf), count, percent);
  }
}

void SimplexMetrics::IncrementAdapterCount(const AdapterMetric& am, const u64 count) { _adapter_counts[am] += count; }

u64 SimplexMetrics::GetAdapterCount(const AdapterMetric& am) const {
  const auto it = _adapter_counts.find(am);
  return it == _adapter_counts.end() ? 0 : it->second;
}

u64 SimplexMetrics::GetSidCollisionCount(const ReadEnd end, const u32 expected_sid, const u32 collided_sid) const {
  const auto it = _sid_collisions.find(SidCollision{end, expected_sid, collided_sid});
  return it == _sid_collisions.end() ? 0 : it->second;
}

f32 SimplexMetrics::CalculatePercentage(const u64 numerator, const u64 denominator) {
  return denominator == 0 ? NAN : static_cast<f32>(static_cast<f64>(numerator) / static_cast<f64>(denominator) * 100.0);
}
}  // namespace xoos::demux
