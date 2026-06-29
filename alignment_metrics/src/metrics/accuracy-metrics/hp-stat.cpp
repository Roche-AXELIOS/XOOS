#include "metrics/accuracy-metrics/hp-stat.h"

#include <filesystem>
#include <fstream>
#include <ranges>
#include <string>
#include <tuple>

#include <fmt/format.h>

#include <csv.hpp>

#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/util/math.h>

#include "alignment-metrics-options.h"
#include "metrics/metrics-names.h"
#include "util/format-util.h"

namespace xoos::alignment_metrics {

const std::string kAllLabel = "any";

HpStats::HpStats(const u8 max_hp_length) {
  for (u8 hp_length = kHpMinLength; hp_length <= max_hp_length; ++hp_length) {
    for (const char base : {'A', 'C', 'G', 'T'}) {
      hp_stats.try_emplace(HpStatKey{base, hp_length});
    }
  }
}

std::vector<std::string> HpStats::GetHeaders() {
  return {kNameHpBase,
          kNameHpLength,
          kNameHpCount,
          kNameHpTotalReads,
          kNameHpSpanningReads,
          kNameHpPercentageSpanning,
          kNameHpMeanSpanningCoverage,
          kNameHpDiscordantReads,
          kNameHpPercentageDiscordant,
          kNameHpLowQualityReads,
          kNameHpPercentageLowQuality,
          kNameHpEffectiveReads,
          kNameHpPercentageEffective,
          kNameHpMeanEffectiveCoverage,
          kNameHpReadsWithInsertion,
          kNameHpReadsWithDeletion,
          kNameHpReadsWithSubstitution,
          kNameHpInsertionErrorRate,
          kNameHpDeletionErrorRate,
          kNameHpIndelErrorRate,
          kNamePhred};
}

// Sums all the hp stats for a given hp_length across all bases
HpStat HpStats::GetTotalHpStat(const u8 hp_length) const {
  HpStat total_hp_stat;
  for (const auto& [hp_stat_key, hp_stat] : hp_stats) {
    if (hp_stat_key.hp_length == hp_length) {
      total_hp_stat.total_reads += hp_stat.total_reads;
      total_hp_stat.spanning_reads += hp_stat.spanning_reads;
      total_hp_stat.discordant_reads += hp_stat.discordant_reads;
      total_hp_stat.low_quality_reads += hp_stat.low_quality_reads;
      total_hp_stat.effective_reads += hp_stat.effective_reads;

      total_hp_stat.count_del += hp_stat.count_del;
      total_hp_stat.count_ins += hp_stat.count_ins;
      total_hp_stat.count_sub += hp_stat.count_sub;
      total_hp_stat.count_indel += hp_stat.count_indel;

      total_hp_stat.hp_count += hp_stat.hp_count;
    }
  }

  return total_hp_stat;
}

void HpStats::IncrementHpCount(const HpStatKey& key) {
  if (hp_stats.contains(key)) {
    ++hp_stats[key].hp_count;
  }
}

void HpStats::WriteTsv(const fs::path& output_path, const io::Comments& comments) const {
  std::ofstream ofstream(output_path);
  auto writer = csv::make_tsv_writer_buffered(ofstream);
  io::WriteTsvComments(writer, comments);
  writer << GetHeaders();
  if (hp_stats.empty()) {
    return;
  }

  // Iterate by hp_length then by base within each length writing out the total for each hp_length
  // after all the bases for that length have been written
  const auto lengths_view = std::ranges::views::keys(hp_stats) |
                            std::ranges::views::transform([](const HpStatKey& key) { return key.hp_length; });
  std::set<u8> unique_lengths(lengths_view.begin(), lengths_view.end());
  for (const u8 hp_length : unique_lengths) {
    for (const char base : {'A', 'C', 'G', 'T'}) {
      HpStatKey key{base, hp_length};
      if (hp_stats.contains(key)) {
        const auto& hp_stat = hp_stats.at(key);
        const f64 mean_spanning_coverage =
            hp_stat.hp_count == 0 ? 0 : static_cast<f64>(hp_stat.spanning_reads) / static_cast<f64>(hp_stat.hp_count);
        const f64 mean_effective_coverage =
            hp_stat.hp_count == 0 ? 0 : static_cast<f64>(hp_stat.effective_reads) / static_cast<f64>(hp_stat.hp_count);
        writer << std::make_tuple(
            std::string(1, base),
            hp_length,
            hp_stat.hp_count,
            hp_stat.total_reads,
            hp_stat.spanning_reads,
            ToPercentageWithPrecision(hp_stat.spanning_reads, hp_stat.total_reads),
            hp_stat.hp_count == 0 ? kNotApplicable : ToStringWithPrecision(mean_spanning_coverage, 0),
            hp_stat.discordant_reads,
            ToPercentageWithPrecision(hp_stat.discordant_reads, hp_stat.total_reads),
            hp_stat.low_quality_reads,
            ToPercentageWithPrecision(hp_stat.low_quality_reads, hp_stat.total_reads),
            hp_stat.effective_reads,
            ToPercentageWithPrecision(hp_stat.effective_reads, hp_stat.total_reads),
            hp_stat.hp_count == 0 ? kNotApplicable : ToStringWithPrecision(mean_effective_coverage, 0),
            hp_stat.count_ins,
            hp_stat.count_del,
            hp_stat.count_sub,
            ToPercentageWithPrecision(hp_stat.count_ins, hp_stat.effective_reads),
            ToPercentageWithPrecision(hp_stat.count_del, hp_stat.effective_reads),
            ToPercentageWithPrecision(hp_stat.count_indel, hp_stat.effective_reads),
            hp_stat.effective_reads == 0
                ? kNotApplicable
                : ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(hp_stat.count_indel) /
                                                               static_cast<f64>(hp_stat.effective_reads))));
      }
    }
    // Write total for this hp_length across all bases
    const auto total_hp_stat = GetTotalHpStat(hp_length);
    const f64 mean_spanning_coverage = total_hp_stat.hp_count == 0 ? 0
                                                                   : static_cast<f64>(total_hp_stat.spanning_reads) /
                                                                         static_cast<f64>(total_hp_stat.hp_count);
    const f64 mean_effective_coverage = total_hp_stat.hp_count == 0 ? 0
                                                                    : static_cast<f64>(total_hp_stat.effective_reads) /
                                                                          static_cast<f64>(total_hp_stat.hp_count);
    writer << std::make_tuple(
        kAllLabel,
        hp_length,
        total_hp_stat.hp_count,
        total_hp_stat.total_reads,
        total_hp_stat.spanning_reads,
        ToPercentageWithPrecision(total_hp_stat.spanning_reads, total_hp_stat.total_reads),
        total_hp_stat.hp_count == 0 ? kNotApplicable : ToStringWithPrecision(mean_spanning_coverage, 0),
        total_hp_stat.discordant_reads,
        ToPercentageWithPrecision(total_hp_stat.discordant_reads, total_hp_stat.total_reads),
        total_hp_stat.low_quality_reads,
        ToPercentageWithPrecision(total_hp_stat.low_quality_reads, total_hp_stat.total_reads),
        total_hp_stat.effective_reads,
        ToPercentageWithPrecision(total_hp_stat.effective_reads, total_hp_stat.total_reads),
        total_hp_stat.hp_count == 0 ? kNotApplicable : ToStringWithPrecision(mean_effective_coverage, 0),
        total_hp_stat.count_ins,
        total_hp_stat.count_del,
        total_hp_stat.count_sub,
        ToPercentageWithPrecision(total_hp_stat.count_ins, total_hp_stat.effective_reads),
        ToPercentageWithPrecision(total_hp_stat.count_del, total_hp_stat.effective_reads),
        ToPercentageWithPrecision(total_hp_stat.count_indel, total_hp_stat.effective_reads),
        total_hp_stat.effective_reads == 0
            ? kNotApplicable
            : ToStringWithPrecision(math::ErrorRateToPhred(static_cast<f64>(total_hp_stat.count_indel) /
                                                           static_cast<f64>(total_hp_stat.effective_reads))));
  }
}

HpStats operator+(const HpStats& lhs, const HpStats& rhs) {
  HpStats result(lhs);
  result += rhs;
  return result;
}

HpStats& HpStats::operator+=(const HpStats& other) {
  for (const auto& [hp_error_key, hp_error] : other.hp_stats) {
    if (hp_stats.contains(hp_error_key)) {
      auto& this_hp_error = hp_stats.at(hp_error_key);
      // base and hp_length should match
      this_hp_error.total_reads += hp_error.total_reads;
      this_hp_error.spanning_reads += hp_error.spanning_reads;
      this_hp_error.discordant_reads += hp_error.discordant_reads;
      this_hp_error.low_quality_reads += hp_error.low_quality_reads;
      this_hp_error.effective_reads += hp_error.effective_reads;
      this_hp_error.count_ins += hp_error.count_ins;
      this_hp_error.count_del += hp_error.count_del;
      this_hp_error.count_indel += hp_error.count_indel;
      this_hp_error.count_sub += hp_error.count_sub;
      this_hp_error.hp_count += hp_error.hp_count;

      // We merge the position-based error count vectors for insertions, deletions, and substitutions
      for (size_t i = 0; i < hp_error.insertion_by_position.size() && i < this_hp_error.insertion_by_position.size();
           ++i) {
        this_hp_error.insertion_by_position[i] += hp_error.insertion_by_position[i];
      }
      for (size_t i = 0; i < hp_error.deletion_by_position.size() && i < this_hp_error.deletion_by_position.size();
           ++i) {
        this_hp_error.deletion_by_position[i] += hp_error.deletion_by_position[i];
      }
      for (size_t i = 0;
           i < hp_error.substitution_by_position.size() && i < this_hp_error.substitution_by_position.size();
           ++i) {
        this_hp_error.substitution_by_position[i] += hp_error.substitution_by_position[i];
      }
      for (size_t i = 0; i < hp_error.indel_by_position.size() && i < this_hp_error.indel_by_position.size(); ++i) {
        this_hp_error.indel_by_position[i] += hp_error.indel_by_position[i];
      }
    }
  }
  return *this;
}

}  // namespace xoos::alignment_metrics
