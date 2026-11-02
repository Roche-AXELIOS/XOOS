#include "generate-pon.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <queue>
#include <ranges>
#include <set>
#include <string>
#include <string_view>
#include <thread>

#include <taskflow/taskflow.hpp>

#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-writer.h>
#include <xoos/log/logging.h>
#include <xoos/util/string-functions.h>

#include "core/vcf-fields.h"
#include "util/file-util.h"

namespace xoos::svc {

// Expected column names in the feature TSV files
static constexpr std::string_view kColChrom = "chrom";
static constexpr std::string_view kColPos = "pos";
static constexpr std::string_view kColRef = "ref";
static constexpr std::string_view kColAlt = "alt";
static constexpr std::string_view kColDuplexAf = "duplex_af";
static constexpr std::string_view kColMapqMean = "mapq_mean";

// Indices of required columns in the feature TSV header
struct ColumnIndices {
  size_t chrom{};
  size_t pos{};
  size_t ref{};
  size_t alt{};
  size_t duplex_af{};
  size_t mapq_mean{};
};

// Per-variant filter thresholds applied while reading feature files
struct PonFilterThresholds {
  f32 min_duplex_af{};
  f32 max_duplex_af{};
  f32 min_mapq_mean{};
};

// ── ShardedVariantMap implementation ──

void ShardedVariantMap::InsertOrUpdate(const VariantId& key, const f32 duplex_af) {
  const auto shard_idx = std::hash<VariantId>{}(key) % kNumShards;
  auto& shard = shards.at(shard_idx);
  const std::lock_guard lock(shard.mutex);

  const auto [it, inserted] = shard.map.try_emplace(
      key,
      VariantStats{
          .occurrence = 1U, .min_duplex_af = duplex_af, .max_duplex_af = duplex_af, .sum_duplex_af = duplex_af});
  if (!inserted) {
    auto& stats = it->second;
    ++stats.occurrence;
    stats.min_duplex_af = std::min(stats.min_duplex_af, duplex_af);
    stats.max_duplex_af = std::max(stats.max_duplex_af, duplex_af);
    stats.sum_duplex_af += duplex_af;
  }
}

size_t ShardedVariantMap::Size() const {
  size_t total = 0;
  for (const auto& shard : shards) {
    total += shard.map.size();
  }
  return total;
}

std::set<std::string, std::less<>> ShardedVariantMap::UniqueContigs() const {
  std::set<std::string, std::less<>> contigs;
  for (const auto& shard : shards) {
    for (const auto& [key, _] : shard.map) {
      contigs.insert(key.chrom);
    }
  }
  return contigs;
}

namespace {
/**
 * @brief Parse the TSV header line and return column indices for required fields.
 * @param header_line Tab-separated header line from the feature file.
 * @param file_path Path to the feature file (for error messages).
 * @return Column indices for chrom, pos, ref, alt, duplex_af, mapq_mean.
 * @throws xoos::Error if any required column is missing.
 */
ColumnIndices ParseHeader(const std::string& header_line, const fs::path& file_path) {
  const auto columns = string::Split(header_line, "\t");
  ColumnIndices indices{};

  struct ColumnMapping {
    std::string_view name;
    size_t* index;
    bool found;
  };

  // Number of required columns: chrom, pos, ref, alt, duplex_af, mapq_mean.
  constexpr size_t kNumRequiredColumns{6};
  std::array<ColumnMapping, kNumRequiredColumns> mappings = {{
      {.name = kColChrom, .index = &indices.chrom, .found = false},
      {.name = kColPos, .index = &indices.pos, .found = false},
      {.name = kColRef, .index = &indices.ref, .found = false},
      {.name = kColAlt, .index = &indices.alt, .found = false},
      {.name = kColDuplexAf, .index = &indices.duplex_af, .found = false},
      {.name = kColMapqMean, .index = &indices.mapq_mean, .found = false},
  }};

  for (size_t i = 0; i < columns.size(); ++i) {
    for (auto& [name, index, found] : mappings) {
      if (columns[i] == name) {
        *index = i;
        found = true;
        break;
      }
    }
  }

  const bool all_found = std::ranges::all_of(mappings, [](const auto& m) { return m.found; });
  if (!all_found) {
    throw error::Error(
        "Feature file {} is missing required columns. Expected: chrom, pos, ref, alt, duplex_af, "
        "mapq_mean",
        file_path);
  }
  return indices;
}

/**
 * @brief Read a single feature TSV file, filter variants, and insert into the sharded map.
 * @param file_path Path to the feature TSV file produced by compute_bam_features.
 * @param thresholds Per-variant filter thresholds (duplex_af range and min mapq_mean).
 * @param sharded_map Shared map to insert passing variants into.
 * @throws xoos::Error if the file cannot be opened or has a malformed header.
 */
void ProcessFeatureFile(const fs::path& file_path,
                        const PonFilterThresholds& thresholds,
                        ShardedVariantMap& sharded_map) {
  std::ifstream input(file_path);
  if (!input.is_open()) {
    throw error::Error("Cannot open feature file: {}", file_path);
  }

  // Find the header line (skip comment lines starting with #)
  std::string line;
  bool header_found = false;
  while (std::getline(input, line)) {
    if (line.empty() || line.starts_with(io::kTsvCommentLinePrefix)) {
      continue;
    }
    header_found = true;
    break;
  }
  if (!header_found) {
    throw error::Error("Cannot find header in feature file: {}", file_path);
  }

  string::Trim(line);
  const auto indices = ParseHeader(line, file_path);
  const auto min_columns =
      std::max({indices.chrom, indices.pos, indices.ref, indices.alt, indices.duplex_af, indices.mapq_mean}) + 1U;

  u64 total_variants = 0;
  u64 passed_variants = 0;

  while (std::getline(input, line)) {
    if (line.empty() || line.starts_with(io::kTsvCommentLinePrefix)) {
      continue;
    }
    string::Trim(line);
    const auto fields = string::Split(line, "\t");
    if (fields.size() < min_columns) {
      Logging::Warn("Skipping malformed line {} in {} (expected {} columns, got {})",
                    total_variants + 1U,
                    file_path,
                    min_columns,
                    fields.size());
      continue;
    }

    ++total_variants;

    f32 duplex_af{};
    f32 mapq_mean{};
    u64 tsv_pos{};
    try {
      duplex_af = std::stof(fields[indices.duplex_af]);
      mapq_mean = std::stof(fields[indices.mapq_mean]);
      tsv_pos = std::stoull(fields[indices.pos]);
    } catch (const std::exception& e) {
      throw error::Error("Failed to parse line {} in {}: {}", total_variants + 1U, file_path, e.what());
    }

    // Apply per-variant filters
    if (duplex_af < thresholds.min_duplex_af || duplex_af > thresholds.max_duplex_af ||
        mapq_mean < thresholds.min_mapq_mean) {
      continue;
    }

    ++passed_variants;

    // Feature TSV positions are 1-based (VariantInfoSerializer writes vid.pos + 1).
    // Convert to 0-based for htslib's SetPosition().
    if (tsv_pos == 0U) {
      throw error::Error("Invalid position 0 at line {} in {} (expected 1-based)", total_variants + 1U, file_path);
    }
    const VariantId key(fields[indices.chrom], tsv_pos - 1U, fields[indices.ref], fields[indices.alt]);

    sharded_map.InsertOrUpdate(key, duplex_af);
  }

  Logging::Info("Processed {}: {} total variants, {} passed filters",
                file_path.filename().string(),
                total_variants,
                passed_variants);
}

/**
 * @brief Build a VCF header with contig lines and PON INFO field definitions.
 *
 * When `command_line` is provided, a `##RocheCommandLine` provenance line is
 * prepended to the header before the contig and INFO entries.
 *
 * @param contigs Sorted set of contig names observed in the variant data.
 * @param command_line Optional CLI metadata (program name, version, arguments);
 *                     when set, emits a ##RocheCommandLine header line.
 * @return Shared pointer to the constructed VCF header.
 */
io::VcfHeaderPtr BuildPonHeader(const std::set<std::string, std::less<>>& contigs,
                                const std::optional<io::CommandLineInfo>& command_line) {
  auto header = io::VcfHeader::Create();

  if (command_line) {
    header->AddCustomMetaDataLine(io::GetCommandInfo(*command_line));
  }

  for (const auto& contig : contigs) {
    header->AddCustomMetaDataLine(fmt::format("##contig=<ID={}>", contig));
  }

  using enum io::FieldType;
  header->AddInfoLine({.id = kFieldPonOcc,
                       .description = "Number of samples containing this variant",
                       .number = io::kNumberOne,
                       .type = kInteger});
  header->AddInfoLine({.id = kFieldPonOccRatio,
                       .description = "Fraction of samples containing this variant",
                       .number = io::kNumberOne,
                       .type = kFloat});
  header->AddInfoLine({.id = kFieldPonMeanDuplexAf,
                       .description = "Mean duplex_af across samples",
                       .number = io::kNumberOne,
                       .type = kFloat});
  header->AddInfoLine({.id = kFieldPonMinDuplexAf,
                       .description = "Minimum duplex_af across samples",
                       .number = io::kNumberOne,
                       .type = kFloat});
  header->AddInfoLine({.id = kFieldPonMaxDuplexAf,
                       .description = "Maximum duplex_af across samples",
                       .number = io::kNumberOne,
                       .type = kFloat});

  header->Sync();
  return header;
}

/**
 * @brief Write a single variant record to the PON VCF.
 * @param writer VCF writer to emit the record.
 * @param key Variant identifier (chrom, pos, ref, alt).
 * @param stats Aggregated statistics for this variant across samples.
 * @param total_samples Total number of input samples (for occurrence ratio).
 */
void WritePonRecord(io::VcfWriter& writer, const VariantId& key, const VariantStats& stats, const u32 total_samples) {
  const auto occ_ratio = static_cast<f32>(stats.occurrence) / static_cast<f32>(total_samples);
  const auto mean_duplex_af = static_cast<f32>(stats.sum_duplex_af / stats.occurrence);

  const io::VcfIntegerFields integer_info = {{kFieldPonOcc, {static_cast<s32>(stats.occurrence)}}};
  const io::VcfFloatFields float_info = {{kFieldPonOccRatio, {occ_ratio}},
                                         {kFieldPonMeanDuplexAf, {mean_duplex_af}},
                                         {kFieldPonMinDuplexAf, {stats.min_duplex_af}},
                                         {kFieldPonMaxDuplexAf, {stats.max_duplex_af}}};
  const io::VcfStringFields string_info = {};

  const io::TypedVcfFields info_fields{
      .integer_fields = integer_info,
      .float_fields = float_info,
      .string_fields = string_info,
      .field_order = std::vector<std::string>{
          kFieldPonOcc, kFieldPonOccRatio, kFieldPonMeanDuplexAf, kFieldPonMinDuplexAf, kFieldPonMaxDuplexAf}};
  const io::TypedVcfFields empty_format{.integer_fields = {}, .float_fields = {}, .string_fields = {}};

  const auto record = writer.CreateRecord(
      key.chrom, static_cast<s32>(key.pos), ".", {key.ref, key.alt}, std::nullopt, "PASS", info_fields, empty_format);
  writer.WriteRecord(record);
}
}  // namespace

void GeneratePon(const GeneratePonParam& param) {
  const auto& feature_files = param.feature_files;
  const auto requested_threads =
      param.threads == 0U ? std::max(1U, std::thread::hardware_concurrency()) : param.threads;
  const auto threads = std::min(requested_threads, feature_files.size());

  Logging::Info("Generating PON from {} feature files using {} thread(s)", feature_files.size(), threads);
  Logging::Info(
      "Filters: duplex_af in [{}, {}], mapq_mean >= {}", param.min_duplex_af, param.max_duplex_af, param.min_mapq_mean);

  // All threads insert into a single sharded map, keeping memory close to single-threaded usage
  ShardedVariantMap sharded_map;

  const PonFilterThresholds thresholds{
      .min_duplex_af = param.min_duplex_af, .max_duplex_af = param.max_duplex_af, .min_mapq_mean = param.min_mapq_mean};

  if (threads <= 1U) {
    for (const auto& file_path : feature_files) {
      ProcessFeatureFile(file_path, thresholds, sharded_map);
    }
  } else {
    tf::Executor executor(threads);
    tf::Taskflow flow;

    for (const auto& file_path : feature_files) {
      flow.emplace([&thresholds, &sharded_map, &fp = file_path] { ProcessFeatureFile(fp, thresholds, sharded_map); });
    }

    executor.run(flow).get();
  }

  const auto total_variants = sharded_map.Size();
  Logging::Info("Total unique variants after filtering: {}", total_variants);

  // Write PON VCF using k-way merge across shards (no second map needed)
  CreateParentDirectoryIfNotExists(param.output_file);
  const auto header = BuildPonHeader(sharded_map.UniqueContigs(), param.command_line);
  io::VcfWriter writer(param.output_file, header);
  writer.WriteHeader();

  // K-way merge: each shard's std::map is sorted by VariantId. Use a min-heap
  // of iterators to emit variants in globally sorted order.
  using ShardIter = VariantMap::const_iterator;

  struct HeapEntry {
    const VariantId* key;
    const VariantStats* stats;
    size_t shard_idx;

    // Order by the pointed-to VariantId so the min-heap emits variants in sorted order.
    auto operator<=>(const HeapEntry& other) const {
      return *key <=> *other.key;
    }
  };

  std::priority_queue<HeapEntry, vec<HeapEntry>, std::greater<>> heap;
  vec<ShardIter> iterators(kNumShards);
  for (size_t i = 0; i < kNumShards; ++i) {
    iterators[i] = sharded_map.shards.at(i).map.cbegin();
    if (iterators[i] != sharded_map.shards.at(i).map.cend()) {
      heap.emplace(&iterators[i]->first, &iterators[i]->second, i);
    }
  }

  const auto total_samples = static_cast<u32>(feature_files.size());
  size_t written = 0;
  while (!heap.empty()) {
    const auto [key, stats, shard_idx] = heap.top();
    heap.pop();
    WritePonRecord(writer, *key, *stats, total_samples);
    ++written;

    ++iterators[shard_idx];
    if (iterators[shard_idx] != sharded_map.shards.at(shard_idx).map.cend()) {
      heap.emplace(&iterators[shard_idx]->first, &iterators[shard_idx]->second, shard_idx);
    }
  }

  writer.Flush();
  Logging::Info("Wrote PON VCF with {} variants to {}", written, param.output_file.string());
}

}  // namespace xoos::svc
