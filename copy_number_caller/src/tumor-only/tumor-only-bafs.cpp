#include "tumor-only/tumor-only-bafs.h"

#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include <fmt/core.h>
#include <fmt/format.h>
#include <htslib/vcf.h>

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-reader.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/file-functions.h>

#include "cli/copy-number-caller-options.h"
#include "io/copy-number-caller-default-filenames.h"
#include "io/write-baf-values.h"
#include "observations/observations.h"
#include "tumor-only/snp-blocklist.h"
#include "util/eigen-math.h"

namespace xoos::cnc {

namespace {

constexpr auto* kPopulationAfField = "AF";
constexpr u64 kProgressInterval = 5'000'000;

struct FilterStats {
  u64 total_sites = 0;
  u64 skipped_not_snp = 0;
  u64 skipped_multiallelic = 0;
  u64 skipped_missing_ad = 0;
  u64 skipped_not_in_pop = 0;
  u64 skipped_low_depth = 0;
  u64 skipped_baf_range = 0;
  u64 skipped_blocklisted = 0;
  u64 kept = 0;
};

void LogFilterStats(const FilterStats& s, const CopyNumberCallerOptions& options) {
  Logging::Info("Tumor-only BAFs summary:");
  Logging::Info("  Total sites in tumor VCF:       {}", s.total_sites);
  Logging::Info("  Skipped (not SNP):              {}", s.skipped_not_snp);
  Logging::Info("  Skipped (multiallelic):         {}", s.skipped_multiallelic);
  Logging::Info("  Skipped (missing/short AD):     {}", s.skipped_missing_ad);
  Logging::Info("  Skipped (low depth < {}):       {}", options.min_tumor_depth, s.skipped_low_depth);
  Logging::Info("  Skipped (BAF outside [{}, {}]): {}",
                options.baf_filter_options.tumor_sample_min_baf,
                options.baf_filter_options.tumor_sample_max_baf,
                s.skipped_baf_range);
  Logging::Info("  Skipped (not in population):    {}", s.skipped_not_in_pop);
  Logging::Info("  Skipped (SNP blocklisted):      {}", s.skipped_blocklisted);
  Logging::Info("  Sites kept:                     {}", s.kept);
}

/// Validates that records within each chromosome arrive in non-decreasing
/// position order.  Throws on violation.
void CheckPositionSorted(const std::string& chrom, const s64 pos, std::string& prev_chrom, s64& prev_pos) {
  if (chrom == prev_chrom && pos < prev_pos) {
    throw std::runtime_error(
        fmt::format("Tumor VCF is not position-sorted: {}:{} appears after {}:{}", chrom, pos, chrom, prev_pos));
  }
  prev_chrom = chrom;
  prev_pos = pos;
}

void BuildObservations(Observations& obvs, const vec<u64>& starts, const vec<u64>& ends, const vec<f64>& ads) {
  obvs.starts = eigen::FromStdVecU64(starts);
  obvs.ends = eigen::FromStdVecU64(ends);
  obvs.obvs = eigen::FromStdVec(ads);
}

/// State for sequential scanning through the population VCF.
struct PopulationScanner {
  // write-through: the scanner advances this caller-owned VcfReader as it scans; a value member would copy a stateful,
  // non-copyable I/O resource.
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
  io::VcfReader& reader;
  std::string current_chrom;
  io::VcfRecordPtr rec;
  bool chrom_exhausted = true;

  explicit PopulationScanner(io::VcfReader& rdr) : reader(rdr) {
  }

  /// Jump to a new chromosome. Sets chrom_exhausted if the chromosome has no records.
  void JumpToChromosome(const std::string& chrom) {
    current_chrom = chrom;
    chrom_exhausted = !reader.SetRegion(chrom, 0, std::numeric_limits<s64>::max());
    if (!chrom_exhausted) {
      Advance();
    }
  }

  /// Advance to the next record. Sets chrom_exhausted if no more records.
  void Advance() {
    rec = reader.GetNextRegionRecord(BCF_UN_INFO);
    chrom_exhausted = (rec == nullptr);
  }

  /// Advance until the current record is at or past the given position.
  void AdvanceTo(const s64 pos) {
    while (!chrom_exhausted && rec->Position() < pos) {
      Advance();
    }
  }
};

/// Check if the population VCF has a matching record at the given position
/// with the same REF/ALT alleles and AF above the threshold.
bool FindPopulationMatch(PopulationScanner& pop,
                         const s64 tumor_pos,
                         const std::string_view tumor_ref,
                         const std::string_view tumor_alt,
                         const f64 min_af) {
  while (!pop.chrom_exhausted && pop.rec->Position() == tumor_pos) {
    if (pop.rec->Allele(0) == tumor_ref && pop.rec->Allele(1) == tumor_alt &&
        pop.rec->HasInfoFieldNoCheck(kPopulationAfField)) {
      auto af_values = pop.rec->GetInfoFieldNoCheck<float>(kPopulationAfField);
      if (!af_values.empty() && static_cast<f64>(af_values[0]) >= min_af) {
        return true;
      }
    }
    pop.Advance();
  }
  return false;
}

/// Resolve the tumor sample index from the VCF header.
s32 GetTumorSampleIndex(const io::VcfHeaderPtr& header, const std::string& tumor_sample_name) {
  const auto sample_indexes = header->GetSampleIndexes();
  const auto it = sample_indexes.find(tumor_sample_name);
  if (it == sample_indexes.end()) {
    throw std::runtime_error("Tumor sample name '" + tumor_sample_name + "' not found in VCF header.");
  }
  return it->second;
}

/// Load the optional SNP blocklist BED, or return an empty blocklist when none is configured.
SnpBlocklist LoadSnpBlocklist(const CopyNumberCallerOptions& options) {
  if (options.snp_blocklist_bed_fname.has_value()) {
    SnpBlocklist blocklist(options.snp_blocklist_bed_fname.value());
    Logging::Info("Loaded SNP blocklist with {} intervals from {}",
                  blocklist.Size(),
                  options.snp_blocklist_bed_fname.value().string());
    return blocklist;
  }
  return {};
}

/// Ref/alt allele depths for a tumor site that passed all tumor-side filters.
struct TumorSiteAd {
  s32 ref_ad{};
  s32 alt_ad{};
};

/// Apply all tumor-side filters (SNP, biallelic, blocklist, AD presence, depth, BAF range) to a
/// single tumor record. Returns the ref/alt ADs when the site passes, or nullopt when it is
/// filtered out, incrementing the corresponding FilterStats counter for the reason.
std::optional<TumorSiteAd> FilterTumorSite(const io::VcfRecord& rec,
                                           const std::string& tumor_chrom,
                                           const s64 tumor_pos,
                                           const s32 tumor_sample_index,
                                           const SnpBlocklist& blocklist,
                                           const CopyNumberCallerOptions& options,
                                           FilterStats& stats) {
  if (!rec.IsSnp()) {
    ++stats.skipped_not_snp;
    return std::nullopt;
  }
  if (ToUnsigned(rec.NumAlleles()) > kBiallelicAlleleCount) {
    ++stats.skipped_multiallelic;
    return std::nullopt;
  }
  // Exclude sites overlapping the SNP blocklist (VCF/BED positions are both 0-based here).
  if (!blocklist.Empty() && blocklist.Contains(tumor_chrom, tumor_pos)) {
    ++stats.skipped_blocklisted;
    return std::nullopt;
  }
  if (!rec.HasFormatFieldNoCheck("AD")) {
    ++stats.skipped_missing_ad;
    return std::nullopt;
  }
  const auto& ad_vals = rec.GetFormatField<s32>("AD");
  const auto idx = static_cast<size_t>(tumor_sample_index) * kBiallelicAlleleCount;
  if (idx + 1 >= ad_vals.size()) {
    ++stats.skipped_missing_ad;
    return std::nullopt;
  }
  const s32 tumor_ref_ad = ad_vals[idx];
  const s32 tumor_alt_ad = ad_vals[idx + 1];
  if (tumor_ref_ad < 0 || tumor_alt_ad < 0) {
    ++stats.skipped_missing_ad;
    return std::nullopt;
  }
  const s64 total_depth = static_cast<s64>(tumor_ref_ad) + static_cast<s64>(tumor_alt_ad);
  if (total_depth <= 0 || total_depth < static_cast<s64>(options.min_tumor_depth)) {
    ++stats.skipped_low_depth;
    return std::nullopt;
  }
  const f64 baf = static_cast<f64>(tumor_alt_ad) / static_cast<f64>(total_depth);
  if (baf < options.baf_filter_options.tumor_sample_min_baf || baf > options.baf_filter_options.tumor_sample_max_baf) {
    ++stats.skipped_baf_range;
    return std::nullopt;
  }
  return TumorSiteAd{.ref_ad = tumor_ref_ad, .alt_ad = tumor_alt_ad};
}

}  // namespace

RefAltObservations ExtractTumorOnlyBafs(const CopyNumberCallerOptions& options) {
  const io::VcfReader tumor_vcf_reader(options.vcf_fname.value().string());  // NOSONAR — mutated by GetNextRecord()
  io::VcfReader pop_vcf_reader(options.population_vcf_fname.value().string());
  if (!pop_vcf_reader.HasIndex()) {
    throw std::runtime_error(
        fmt::format("Population VCF must be bgzip-compressed and indexed (e.g. .tbi/.csi) for region queries: {}",
                    options.population_vcf_fname.value().string()));
  }
  const auto tumor_sample_index =
      GetTumorSampleIndex(tumor_vcf_reader.GetHeader(), options.sample_metadata_options.tumor_sample_name.value());
  const SnpBlocklist blocklist = LoadSnpBlocklist(options);

  RefAltObservations result;
  Observations& ref_obvs = result.ref_obvs;
  Observations& alt_obvs = result.alt_obvs;
  vec<u64> starts;
  vec<u64> ends;
  vec<f64> ref_ads;
  vec<f64> alt_ads;
  FilterStats stats;
  PopulationScanner pop(pop_vcf_reader);
  std::string prev_chrom;
  s64 prev_pos = -1;
  io::VcfRecordPtr rec;
  while ((rec = tumor_vcf_reader.GetNextRecord()) != nullptr) {
    ++stats.total_sites;

    if (stats.total_sites % kProgressInterval == 0) {
      Logging::Info("  Progress: {} tumor sites processed, {} kept so far", stats.total_sites, stats.kept);
    }

    const auto& tumor_chrom = rec->Chromosome();
    const auto tumor_pos = rec->Position();
    CheckPositionSorted(tumor_chrom, tumor_pos, prev_chrom, prev_pos);

    const auto tumor_ad = FilterTumorSite(*rec, tumor_chrom, tumor_pos, tumor_sample_index, blocklist, options, stats);
    if (!tumor_ad.has_value()) {
      continue;
    }
    const s32 tumor_ref_ad = tumor_ad->ref_ad;
    const s32 tumor_alt_ad = tumor_ad->alt_ad;

    // On chromosome transition, jump the population reader
    if (tumor_chrom != pop.current_chrom) {
      pop.JumpToChromosome(tumor_chrom);
    }

    if (pop.chrom_exhausted) {
      ++stats.skipped_not_in_pop;
      continue;
    }
    pop.AdvanceTo(tumor_pos);
    if (!FindPopulationMatch(pop, tumor_pos, rec->Allele(0), rec->Allele(1), options.min_population_af)) {
      ++stats.skipped_not_in_pop;
      continue;
    }

    // Keep this site
    ref_obvs.contigs.push_back(tumor_chrom);
    alt_obvs.contigs.push_back(tumor_chrom);
    starts.push_back(static_cast<u64>(tumor_pos));
    ends.push_back(static_cast<u64>(tumor_pos + 1));
    ref_ads.push_back(static_cast<f64>(tumor_ref_ad));
    alt_ads.push_back(static_cast<f64>(tumor_alt_ad));
    ++stats.kept;
  }

  BuildObservations(ref_obvs, starts, ends, ref_ads);
  BuildObservations(alt_obvs, starts, ends, alt_ads);

  LogFilterStats(stats, options);
  return result;
}

void TumorOnlyBafsMain(const CopyNumberCallerOptions& options) {
  const auto baf_out = options.output_dir / kDefaultBafOutput;
  file::CheckOutputFilePathsForWriting({baf_out}, options.overwrite);

  const RefAltObservations obvs = ExtractTumorOnlyBafs(options);

  Logging::Info("Writing BAF values to {}", baf_out.string());
  WriteBafValues(obvs.ref_obvs, obvs.alt_obvs, baf_out, options.command_line_info);
}

}  // namespace xoos::cnc
