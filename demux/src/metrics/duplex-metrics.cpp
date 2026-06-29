#include "metrics/duplex-metrics.h"

#include <fmt/format.h>
#include <xoos/histogram/histogram-summary.h>
#include <xoos/log/logging.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <functional>
#include <numeric>
#include <string>

#include "core/demux-and-trim-pipeline.h"
#include "csv.hpp"
#include "metric-filenames.h"
#include "metrics-constraints.h"
#include "metrics.h"
#include "xoos/types/str-container.h"

namespace xoos::demux {

thread_local concurrent::EnumerableThreadLocal<DuplexMetrics> DuplexMetrics::instance{
    std::make_shared<DuplexMetrics>()};

DuplexMetrics::DuplexMetrics()
    : total_length_distr(metrics_constraints::max_logged_read_length + 1, 0),
      unassigned_length_distr(metrics_constraints::max_logged_read_length + 1, 0),
      no_hairpin_length_distr(metrics_constraints::max_logged_read_length + 1, 0),
      full_duplex_length_distr(metrics_constraints::max_sid_id_index + 1,
                               LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)),
      partial_duplex_length_distr(metrics_constraints::max_sid_id_index + 1,
                                  LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)),
      endadapter_position_distr(metrics_constraints::max_sid_id_index + 1,
                                LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)),
      passing_length_distr(metrics_constraints::max_sid_id_index + 1,
                           LengthHistogram(metrics_constraints::max_logged_read_length + 1, 0)) {}

static void AddOther(PerSIDHistogram& hist, const PerSIDHistogram& other) {
  for (size_t i = 0; i < hist.size(); ++i) {
    hist[i] += other[i];
  }
}

static void AddOther(PerSIDCount& count, const PerSIDCount& other) {
  for (size_t i = 0; i < count.size(); ++i) {
    count[i] += other[i];
  }
}

DuplexMetrics& DuplexMetrics::Instance() { return *instance.Local(); }

void DuplexMetrics::Add(const DuplexMetrics& other) {
  // merge unassigned counts
  unassigned_counts.read_too_long += other.unassigned_counts.read_too_long;
  unassigned_counts.no_hairpin_found += other.unassigned_counts.no_hairpin_found;
  unassigned_counts.read_too_short += other.unassigned_counts.read_too_short;

  // merge failed assigned counts
  AddOther(failed_assigned_counts.consensus_too_long, other.failed_assigned_counts.consensus_too_long);
  AddOther(failed_assigned_counts.too_many_errors, other.failed_assigned_counts.too_many_errors);
  AddOther(failed_assigned_counts.trimmed_read_too_short, other.failed_assigned_counts.trimmed_read_too_short);
  AddOther(failed_assigned_counts.failed_hairpin_stem_trim_reads,
           other.failed_assigned_counts.failed_hairpin_stem_trim_reads);

  // merge mid adapter counts
  AddOther(midadapter_counts.found_by_global_symmetry, other.midadapter_counts.found_by_global_symmetry);
  AddOther(midadapter_counts.found_by_local_symmetry, other.midadapter_counts.found_by_local_symmetry);
  AddOther(midadapter_counts.found_by_string_compare, other.midadapter_counts.found_by_string_compare);

  // merge strand counts
  AddOther(strand_counts.fw, other.strand_counts.fw);
  AddOther(strand_counts.rv, other.strand_counts.rv);
  AddOther(strand_counts.fw_sig, other.strand_counts.fw_sig);
  AddOther(strand_counts.rv_sig, other.strand_counts.rv_sig);

  // merge passing read counts
  AddOther(passing_counts.full_duplex, other.passing_counts.full_duplex);
  AddOther(passing_counts.partial_duplex, other.passing_counts.partial_duplex);
  AddOther(passing_counts.longer_r2, other.passing_counts.longer_r2);
  AddOther(passing_counts.longer_r2_full_duplex, other.passing_counts.longer_r2_full_duplex);
  AddOther(passing_counts.both_umi, other.passing_counts.both_umi);
  AddOther(passing_counts.only_5p_umi, other.passing_counts.only_5p_umi);
  AddOther(passing_counts.only_3p_umi, other.passing_counts.only_3p_umi);
  AddOther(passing_counts.no_endadapter, other.passing_counts.no_endadapter);

  // merge base counts
  AddOther(base_counts.concordant, other.base_counts.concordant);
  AddOther(base_counts.discordant, other.base_counts.discordant);
  AddOther(base_counts.simplex, other.base_counts.simplex);

  // merge read length distributions
  AddOther(full_duplex_length_distr, other.full_duplex_length_distr);
  AddOther(partial_duplex_length_distr, other.partial_duplex_length_distr);
  AddOther(endadapter_position_distr, other.endadapter_position_distr);
  AddOther(passing_length_distr, other.passing_length_distr);
  unassigned_length_distr += other.unassigned_length_distr;
  no_hairpin_length_distr += other.no_hairpin_length_distr;
  total_length_distr += other.total_length_distr;
}

// For early stopping, returns the count of the sample with the smallest count of concordant duplex bases
size_t DuplexMetrics::MinConcordDupBases() {
  PerSIDCount total(metrics_constraints::max_sid_id_index + 1, 0);
  // Sum all the concordant_duplex_bases_total per sample
  const std::function func = [&total](const DuplexMetrics& item) { AddOther(total, item.base_counts.concordant); };
  instance.ForEach(func);
  // Return the count of the sample with the minimum concordant duplex bases
  return *std::ranges::min_element(total.begin(), total.end());
}

DuplexMetrics DuplexMetrics::SumTotal() {
  DuplexMetrics total;
  const std::function func = [&total](const DuplexMetrics& item) { total.Add(item); };
  instance.ForEach(func);
  return total;
}

/**
 * Maps SID names to their corresponding IDs using the provided SID pool. This is useful for generating human-readable
 * metrics outputs.
 * @param sid_pool  The SID pool containing mappings of SID IDs to their corresponding names and other metadata.
 * @return  A map where the keys are SID names and the values are their corresponding IDs. If a SID name is not found in
 * the pool, it will not be included in the map.
 */
StrMap<u32> ConvertSidPoolToSidNameToIdMap(const detail::SidPool& sid_pool) {
  StrMap<u32> sid_name_to_id;
  for (const auto& [_, sid_id] : sid_pool) {
    sid_name_to_id[sid_id.name] = sid_id.id;
  }
  return sid_name_to_id;
}

void DuplexMetrics::WriteMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool) const {
  fs::create_directory(params.out_dir / kMetricsDirectory);

  const StrMap<u32> sid_name_to_id = ConvertSidPoolToSidNameToIdMap(sid_pool);
  WriteRunMetrics(params, sid_pool, sid_name_to_id);
  WriteSampleMetrics(params, sid_pool, sid_name_to_id);
  WriteReadLengthDistributions(params);
  WriteSampleAssignmentMetrics(params, sid_pool);
}

static float Percentage(u64 count, u64 div) {
  float percentage = 0.0f;
  if (div > 0) {
    percentage = (static_cast<float>(count) / static_cast<float>(div)) * 100.0f;
  }
  return percentage;
}

static f64 ComputeMeanReadLength(const PerSIDHistogram& passing_length_distr) {
  u64 total_count = 0;
  u64 total_read_length_sum = 0;
  for (auto& sid_hist : passing_length_distr) {
    total_count += ComputeCount(sid_hist);
    total_read_length_sum += histogram::ComputeValueSum(sid_hist);
  }
  return total_count > 0 ? static_cast<f64>(total_read_length_sum) / static_cast<f64>(total_count) : 0.0f;
}

std::unordered_set<u32> DuplexMetrics::FindAssignedSids(const StrMap<u32>& sid_name_to_id) const {
  // Find all SIDs that have at least one read assigned to them
  std::unordered_set<u32> found_sids;
  for (const auto& [name, id] : sid_name_to_id) {
    if (id < passing_counts.full_duplex.size() &&
        passing_counts.full_duplex.at(id) + passing_counts.partial_duplex.at(id) +
                failed_assigned_counts.consensus_too_long.at(id) + failed_assigned_counts.too_many_errors.at(id) +
                failed_assigned_counts.trimmed_read_too_short.at(id) >
            0) {
      found_sids.insert(id);
    }
  }
  return found_sids;
}

void DuplexMetrics::WriteRunMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool,
                                    const StrMap<u32>& sid_name_to_id) const {
  // compute composite metrics
  // total assigned reads = passing reads + failed assigned reads
  // failed assigned reads metrics
  const auto consensus_too_long = accumulate(failed_assigned_counts.consensus_too_long.begin(),
                                             failed_assigned_counts.consensus_too_long.end(), 0UL);
  const auto too_many_errors =
      accumulate(failed_assigned_counts.too_many_errors.begin(), failed_assigned_counts.too_many_errors.end(), 0UL);
  const auto trimmed_read_too_short = accumulate(failed_assigned_counts.trimmed_read_too_short.begin(),
                                                 failed_assigned_counts.trimmed_read_too_short.end(), 0UL);
  const auto failed_midadapter_trim_reads{accumulate(failed_assigned_counts.failed_hairpin_stem_trim_reads.begin(),
                                                     failed_assigned_counts.failed_hairpin_stem_trim_reads.end(), 0UL)};
  const auto failed_assigned =
      consensus_too_long + too_many_errors + trimmed_read_too_short + failed_midadapter_trim_reads;

  // passing reads = full duplex reads + partial duplex reads
  const auto full_duplex_reads =
      std::accumulate(passing_counts.full_duplex.begin(), passing_counts.full_duplex.end(), 0UL);
  const auto partial_duplex_reads =
      std::accumulate(passing_counts.partial_duplex.begin(), passing_counts.partial_duplex.end(), 0UL);
  const auto passing_reads = full_duplex_reads + partial_duplex_reads;

  const auto unassigned_reads =
      unassigned_counts.read_too_short + unassigned_counts.read_too_long + unassigned_counts.no_hairpin_found;

  // total assigned reads = passing reads + failed assigned reads
  const auto assigned_reads = failed_assigned + passing_reads;
  const auto total_reads = assigned_reads + unassigned_reads;

  // mid adapter finding counts
  const auto found_string(std::accumulate(midadapter_counts.found_by_string_compare.begin(),
                                          midadapter_counts.found_by_string_compare.end(), 0UL));
  const auto found_global(std::accumulate(midadapter_counts.found_by_global_symmetry.begin(),
                                          midadapter_counts.found_by_global_symmetry.end(), 0UL));
  const auto found_local(std::accumulate(midadapter_counts.found_by_local_symmetry.begin(),
                                         midadapter_counts.found_by_local_symmetry.end(), 0UL));

  // compute passing read counts
  const auto longer_r2{std::accumulate(passing_counts.longer_r2.begin(), passing_counts.longer_r2.end(), 0UL)};
  const auto longer_r2_full_duplex{
      std::accumulate(passing_counts.longer_r2_full_duplex.begin(), passing_counts.longer_r2_full_duplex.end(), 0UL)};
  const auto both_umi{std::accumulate(passing_counts.both_umi.begin(), passing_counts.both_umi.end(), 0UL)};
  const auto only_5p_umi{std::accumulate(passing_counts.only_5p_umi.begin(), passing_counts.only_5p_umi.end(), 0UL)};
  const auto only_3p_umi{std::accumulate(passing_counts.only_3p_umi.begin(), passing_counts.only_3p_umi.end(), 0UL)};
  const auto no_endadapter{accumulate(passing_counts.no_endadapter.begin(), passing_counts.no_endadapter.end(), 0UL)};

  // compute strand counts
  const auto strand_fw_reads{std::accumulate(strand_counts.fw.begin(), strand_counts.fw.end(), 0UL)};
  const auto strand_rv_reads{std::accumulate(strand_counts.rv.begin(), strand_counts.rv.end(), 0UL)};
  const auto strand_fw_sig_reads{std::accumulate(strand_counts.fw_sig.begin(), strand_counts.fw_sig.end(), 0UL)};
  const auto strand_rv_sig_reads{std::accumulate(strand_counts.rv_sig.begin(), strand_counts.rv_sig.end(), 0UL)};

  // compute base counts
  // concordant bases = concordant duplex bases + discordant bases + non-duplex bases
  const auto concordant{std::accumulate(base_counts.concordant.begin(), base_counts.concordant.end(), 0UL)};
  const auto discordant{std::accumulate(base_counts.discordant.begin(), base_counts.discordant.end(), 0UL)};
  const auto nonduplex{std::accumulate(base_counts.simplex.begin(), base_counts.simplex.end(), 0UL)};
  const auto total_bases = concordant + discordant + nonduplex;

  const auto out_dir = params.out_dir / kMetricsDirectory;
  const auto out_file_name{out_dir / kRunMetricsFile};

  std::ofstream out_file(out_file_name);

  // Set floating point precision for the entire file stream
  out_file << std::fixed << std::setprecision(2);

  // Create a TSV writer from the file stream
  auto writer = csv::make_tsv_writer(out_file);
  io::WriteTsvComments(writer, params.comments);

  // Write the header row
  writer << std::vector<std::string>{"metrics_name", "count", "percentage"};

  // Write data rows using std::make_tuple
  writer << std::make_tuple(kTotalReads, total_reads, 100.00);
  writer << std::make_tuple(kAssignedReads, assigned_reads, Percentage(assigned_reads, total_reads));

  // passing read metrics
  writer << std::make_tuple(kPassingReads, passing_reads, Percentage(passing_reads, total_reads));
  writer << std::make_tuple(kFullDuplexReads, full_duplex_reads, Percentage(full_duplex_reads, total_reads));
  writer << std::make_tuple(kPartialDuplexReads, partial_duplex_reads, Percentage(partial_duplex_reads, total_reads));

  // unassigned read metrics
  writer << std::make_tuple(kUnassignedReads, unassigned_reads, Percentage(unassigned_reads, total_reads));
  writer << std::make_tuple(kNoHairpinReads, unassigned_counts.no_hairpin_found,
                            Percentage(unassigned_counts.no_hairpin_found, total_reads));
  writer << std::make_tuple(kTooLongReads, unassigned_counts.read_too_long,
                            Percentage(unassigned_counts.read_too_long, total_reads));
  writer << std::make_tuple(kTooShortReads, unassigned_counts.read_too_short,
                            Percentage(unassigned_counts.read_too_short, total_reads));

  // failed assigned read metrics
  writer << std::make_tuple(kFailedAssignedReads, failed_assigned, Percentage(failed_assigned, total_reads));
  writer << std::make_tuple(kTooManyErrorsReads, too_many_errors, Percentage(too_many_errors, total_reads));
  writer << std::make_tuple(kTooShortTrimmedReads, trimmed_read_too_short,
                            Percentage(trimmed_read_too_short, total_reads));
  writer << std::make_tuple(kTooLongConsensusReads, consensus_too_long, Percentage(consensus_too_long, total_reads));
  writer << std::make_tuple(kFailedHairpinStemTrimReads, failed_midadapter_trim_reads,
                            Percentage(failed_midadapter_trim_reads, total_reads));

  // mid adapter finding metrics
  writer << std::make_tuple(kFoundByStringCompare, found_string, Percentage(found_string, assigned_reads));
  writer << std::make_tuple(kFoundByGlobalSymmetry, found_global, Percentage(found_global, assigned_reads));
  writer << std::make_tuple(kFoundByLocalSymmetry, found_local, Percentage(found_local, assigned_reads));

  // passing read metrics do not contribute to passing read counts (can co-occur)
  writer << std::make_tuple(kLongerR2Reads, longer_r2, Percentage(longer_r2, passing_reads));
  writer << std::make_tuple(kLongerR2FullDuplexReads, longer_r2_full_duplex,
                            Percentage(longer_r2_full_duplex, passing_reads));
  writer << std::make_tuple(kBothUmiReads, both_umi, Percentage(both_umi, passing_reads));
  writer << std::make_tuple(k5pUmiReads, only_5p_umi, Percentage(only_5p_umi, passing_reads));
  writer << std::make_tuple(k3pUmiReads, only_3p_umi, Percentage(only_3p_umi, passing_reads));
  writer << std::make_tuple(kNoEndadapterReads, no_endadapter, Percentage(no_endadapter, passing_reads));

  // strand metrics
  writer << std::make_tuple(kStrandFwReads, strand_fw_reads, Percentage(strand_fw_reads, passing_reads));
  writer << std::make_tuple(kStrandRvReads, strand_rv_reads, Percentage(strand_rv_reads, passing_reads));
  writer << std::make_tuple(kStrandFwSigReads, strand_fw_sig_reads, Percentage(strand_fw_sig_reads, passing_reads));
  writer << std::make_tuple(kStrandRvSigReads, strand_rv_sig_reads, Percentage(strand_rv_sig_reads, passing_reads));

  // base counts
  writer << std::make_tuple(kTotalBases, total_bases, 100.00);
  writer << std::make_tuple(kConcordantDuplexBases, concordant, Percentage(concordant, total_bases));
  writer << std::make_tuple(kDiscordantDuplexBases, discordant, Percentage(discordant, total_bases));
  writer << std::make_tuple(kDuplexBases, concordant + discordant, Percentage(concordant + discordant, total_bases));
  writer << std::make_tuple(kNonDuplexBases, nonduplex, Percentage(nonduplex, total_bases));

  // sids
  writer << std::make_tuple(kNumExpectedSids, sid_pool.size(), 100.00);

  const auto num_sids = FindAssignedSids(sid_name_to_id);
  writer << std::make_tuple(kNumSids, num_sids.size(), Percentage(num_sids.size(), sid_pool.size()));

  auto mean_read_length = ComputeMeanReadLength(passing_length_distr);
  std::string mean_read_length_str = fmt::format("{:.0f}", mean_read_length);
  writer << std::make_tuple(kMeanPassingReadLength, mean_read_length_str, "100.00");
  out_file.close();
}

static std::vector<std::string> ToStringVector(const StrMap<u32>& sid_name_to_id, const std::string_view metric,
                                               const PerSIDCount& count) {
  std::vector<std::string> result;
  result.reserve(sid_name_to_id.size() + 1);
  result.emplace_back(metric.data());
  for (const auto& [name, id] : sid_name_to_id) {
    if (id < count.size()) {
      result.emplace_back(std::to_string(count.at(id)));
    } else {
      result.emplace_back("0");
    }
  }
  return result;
}

template <typename ComputeFunc>
static std::vector<std::string> ComputePerSidMetric(const StrMap<u32>& sid_name_to_id, const std::string_view metric,
                                                    const ComputeFunc& compute_value) {
  std::vector<std::string> result;
  result.reserve(sid_name_to_id.size() + 1);
  result.emplace_back(metric.data());
  for (const auto& [name, id] : sid_name_to_id) {
    const u64 value = compute_value(id);
    result.emplace_back(std::to_string(value));
  }
  return result;
}

void DuplexMetrics::WriteSampleMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool,
                                       const StrMap<u32>& sid_name_to_id) const {
  const auto out_file_name{params.out_dir / kMetricsDirectory / kSampleMetricsFile};

  std::ofstream out_file(out_file_name);
  // Create a TSV writer from the file stream
  auto writer = csv::make_tsv_writer(out_file);
  io::WriteTsvComments(writer, params.comments);

  // write the header row
  std::vector<std::string> current_row{kMetric};
  current_row.reserve(sid_name_to_id.size() + 1);
  for (const auto& [name, id] : sid_name_to_id) {
    current_row.emplace_back(name);
  }
  writer << current_row;

  // index sequence row
  current_row.assign(1, kIndexSequence);
  for (const auto& [name, id] : sid_name_to_id) {
    current_row.emplace_back(sid_pool.at(id).sequence);
  }
  writer << current_row;

  // assigned reads
  writer << ComputePerSidMetric(sid_name_to_id, kAssignedReads, [this](const u32 id) -> u64 {
    if (id < passing_counts.full_duplex.size()) {
      return passing_counts.full_duplex.at(id) + passing_counts.partial_duplex.at(id) +
             failed_assigned_counts.consensus_too_long.at(id) + failed_assigned_counts.too_many_errors.at(id) +
             failed_assigned_counts.trimmed_read_too_short.at(id) +
             failed_assigned_counts.failed_hairpin_stem_trim_reads.at(id);
    }
    return 0;
  });

  // passing reads and related counts
  writer << ComputePerSidMetric(sid_name_to_id, kPassingReads, [this](const u32 id) -> u64 {
    if (id < passing_counts.full_duplex.size()) {
      return passing_counts.full_duplex.at(id) + passing_counts.partial_duplex.at(id);
    }
    return 0;
  });

  writer << ToStringVector(sid_name_to_id, kFullDuplexReads, passing_counts.full_duplex);
  writer << ToStringVector(sid_name_to_id, kPartialDuplexReads, passing_counts.partial_duplex);

  // failed assigned reads
  writer << ComputePerSidMetric(sid_name_to_id, kFailedAssignedReads, [this](const u32 id) -> u64 {
    if (id < failed_assigned_counts.consensus_too_long.size()) {
      return failed_assigned_counts.consensus_too_long.at(id) + failed_assigned_counts.too_many_errors.at(id) +
             failed_assigned_counts.trimmed_read_too_short.at(id) +
             failed_assigned_counts.failed_hairpin_stem_trim_reads.at(id);
    }
    return 0;
  });

  // failed assigned read counts
  writer << ToStringVector(sid_name_to_id, kTooManyErrorsReads, failed_assigned_counts.too_many_errors);
  writer << ToStringVector(sid_name_to_id, kTooShortTrimmedReads, failed_assigned_counts.trimmed_read_too_short);
  writer << ToStringVector(sid_name_to_id, kTooLongConsensusReads, failed_assigned_counts.consensus_too_long);
  writer << ToStringVector(sid_name_to_id, kFailedHairpinStemTrimReads,
                           failed_assigned_counts.failed_hairpin_stem_trim_reads);

  // midadapter finding method counts
  writer << ToStringVector(sid_name_to_id, kFoundByStringCompare, midadapter_counts.found_by_string_compare);
  writer << ToStringVector(sid_name_to_id, kFoundByGlobalSymmetry, midadapter_counts.found_by_global_symmetry);
  writer << ToStringVector(sid_name_to_id, kFoundByLocalSymmetry, midadapter_counts.found_by_local_symmetry);

  // passing reads with unique properties
  writer << ToStringVector(sid_name_to_id, kLongerR2Reads, passing_counts.longer_r2);
  writer << ToStringVector(sid_name_to_id, kLongerR2FullDuplexReads, passing_counts.longer_r2_full_duplex);
  // Add these lines where other passing read metrics are written
  writer << ToStringVector(sid_name_to_id, kBothUmiReads, passing_counts.both_umi);
  writer << ToStringVector(sid_name_to_id, k5pUmiReads, passing_counts.only_5p_umi);
  writer << ToStringVector(sid_name_to_id, k3pUmiReads, passing_counts.only_3p_umi);
  writer << ToStringVector(sid_name_to_id, kNoEndadapterReads, passing_counts.no_endadapter);

  // strand metrics
  writer << ToStringVector(sid_name_to_id, kStrandFwReads, strand_counts.fw);
  writer << ToStringVector(sid_name_to_id, kStrandRvReads, strand_counts.rv);
  writer << ToStringVector(sid_name_to_id, kStrandFwSigReads, strand_counts.fw_sig);
  writer << ToStringVector(sid_name_to_id, kStrandRvSigReads, strand_counts.rv_sig);

  // base counts
  writer << ComputePerSidMetric(sid_name_to_id, kTotalBases, [this](const u32 id) -> u64 {
    if (id < base_counts.concordant.size()) {
      return base_counts.concordant.at(id) + base_counts.discordant.at(id) + base_counts.simplex.at(id);
    }
    return 0;
  });

  writer << ToStringVector(sid_name_to_id, kConcordantDuplexBases, base_counts.concordant);
  writer << ToStringVector(sid_name_to_id, kDiscordantDuplexBases, base_counts.discordant);

  // create a PerSid Count for duplex bases by summing concordant and discordant counts
  writer << ComputePerSidMetric(sid_name_to_id, kDuplexBases, [this](const u32 id) -> u64 {
    if (id < base_counts.concordant.size() && id < base_counts.discordant.size()) {
      return base_counts.concordant.at(id) + base_counts.discordant.at(id);
    }
    return 0;
  });

  writer << ToStringVector(sid_name_to_id, kNonDuplexBases, base_counts.simplex);

  writer << ComputePerSidMetric(sid_name_to_id, kMeanPassingReadLength, [this](const u32 id) -> u64 {
    if (id < passing_length_distr.size()) {
      const auto histo = passing_length_distr.at(id);
      return histo.IsEmpty() ? 0 : ComputeMean(passing_length_distr.at(id));
    }
    return 0;
  });

  out_file.close();
}

static void WriteLengthHistogram(const LengthHistogram& histo, const std::string& title, const fs::path& out_file_name,
                                 const io::Comments& comments) {
  histogram::Histograms<u64> histograms;
  histograms.emplace_back(title, histo);
  histogram::WriteHistograms(histograms, out_file_name, "length", histogram::kMaxLastBin, {}, comments);
}

static void WriteSampleLengthHistogram(const PerSIDHistogram& histo, const DuplexMetrics::SidPool& sid_pool,
                                       const fs::path& out_file_name, const io::Comments& comments) {
  histogram::Histograms<u64> histograms;
  histograms.reserve(sid_pool.size());
  for (const auto& [_, sid] : sid_pool) {
    histograms.emplace_back(sid.name, histo[sid.id]);
  }
  histogram::WriteHistograms(histograms, out_file_name, "length", histogram::kMaxLastBin, {}, comments);
}

void DuplexMetrics::WriteReadLengthDistributions(const DemuxAndTrimParam& params) const {
  const auto out_dir = params.out_dir / kMetricsDirectory;
  const auto& comments = params.comments;
  WriteLengthHistogram(total_length_distr, "total", out_dir / kTotalReadLengthDistr, comments);
  WriteLengthHistogram(unassigned_length_distr, "unassigned", out_dir / kUnassignedReadLengthDistr, comments);
  WriteLengthHistogram(no_hairpin_length_distr, "no_hairpin", out_dir / kNoHairpinReadLengthDistr, comments);
}

void DuplexMetrics::WriteSampleAssignmentMetrics(const DemuxAndTrimParam& params, const SidPool& sid_pool) const {
  const auto out_dir = params.out_dir / kMetricsDirectory;
  const auto& comments = params.comments;
  WriteSampleLengthHistogram(full_duplex_length_distr, sid_pool, out_dir / kFullDuplexReadLengthDistr, comments);
  WriteSampleLengthHistogram(partial_duplex_length_distr, sid_pool, out_dir / kPartialDuplexReadLengthDistr, comments);
  WriteSampleLengthHistogram(endadapter_position_distr, sid_pool, out_dir / kEndadapterPositionDistr, comments);
  WriteSampleLengthHistogram(passing_length_distr, sid_pool, out_dir / kPassingReadLengthDistr, comments);
}

void ReportStrandMetrics(const DuplexMetrics& global_results) {
  const auto strand_fw_reads{
      std::accumulate(global_results.strand_counts.fw.begin(), global_results.strand_counts.fw.end(), 0UL)};
  const auto strand_rv_reads{
      std::accumulate(global_results.strand_counts.rv.begin(), global_results.strand_counts.rv.end(), 0UL)};
  const auto strand_fw_sig_reads{
      std::accumulate(global_results.strand_counts.fw_sig.begin(), global_results.strand_counts.fw_sig.end(), 0UL)};
  const auto strand_rv_sig_reads{
      std::accumulate(global_results.strand_counts.rv_sig.begin(), global_results.strand_counts.rv_sig.end(), 0UL)};
  const auto total_fw_reads = strand_fw_reads + strand_fw_sig_reads;
  const auto total_rv_reads = strand_rv_reads + strand_rv_sig_reads;
  const auto full_duplex_reads = std::accumulate(global_results.passing_counts.full_duplex.begin(),
                                                 global_results.passing_counts.full_duplex.end(), 0UL);
  const auto partial_duplex_reads = std::accumulate(global_results.passing_counts.partial_duplex.begin(),
                                                    global_results.passing_counts.partial_duplex.end(), 0UL);
  const auto passing_reads = full_duplex_reads + partial_duplex_reads;
  const auto ambiguous_reads = passing_reads - total_fw_reads - total_rv_reads;
  const auto fw_percentage = static_cast<double>(total_fw_reads) / static_cast<double>(passing_reads);
  const auto rv_percentage = static_cast<double>(total_rv_reads) / static_cast<double>(passing_reads);
  const auto ambiguous_reads_percentage = static_cast<double>(ambiguous_reads) / static_cast<double>(passing_reads);
  Logging::Info(
      "Strand detection: Forward reads: {} ({:.2f}%), Reverse reads: {} ({:.2f}%), "
      "Ambiguous reads: {} ({:.2f}%)",
      total_fw_reads, fw_percentage * 100.0, total_rv_reads, rv_percentage * 100.0, ambiguous_reads,
      ambiguous_reads_percentage * 100.0);
}

}  // namespace xoos::demux
