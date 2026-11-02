#pragma once

#include <array>
#include <functional>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <string>

#include <xoos/io/metadata-util.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/variant-id.h"

namespace xoos::svc {

/** @brief Accumulated statistics for a single variant across all samples. */
struct VariantStats {
  u32 occurrence{0};
  f32 min_duplex_af{};
  f32 max_duplex_af{};
  f64 sum_duplex_af{};
};

using VariantMap = std::map<VariantId, VariantStats>;

// Number of shards for the concurrent variant map. Chosen to be large enough
// to minimize lock contention across threads while keeping overhead low.
static constexpr size_t kNumShards = 128;

/** @brief A single shard: a mutex protecting its own VariantMap partition. */
struct VariantShard {
  std::mutex mutex;
  VariantMap map;
};

/**
 * @brief Thread-safe sharded variant map.
 * @details Variants are distributed across shards by hashing the variant key,
 *          so concurrent threads rarely contend on the same lock.
 */
struct ShardedVariantMap {
  std::array<VariantShard, kNumShards> shards;

  /**
   * @brief Insert or update a variant in the appropriate shard.
   * @param key Variant identifier (chrom, pos, ref, alt).
   * @param duplex_af Duplex allele fraction for this observation.
   * @details Thread-safe. Selects a shard via hash of the key.
   */
  void InsertOrUpdate(const VariantId& key, f32 duplex_af);

  /**
   * @brief Total number of unique variants across all shards.
   * @return Sum of sizes of all shard maps.
   */
  size_t Size() const;

  /**
   * @brief Collect unique chromosome names across all shards.
   * @return Sorted set of contig names.
   */
  std::set<std::string, std::less<>> UniqueContigs() const;
};

/** @brief CLI parameters for PON generation. */
struct GeneratePonParam {
  vec<fs::path> feature_files;
  fs::path output_file;
  size_t threads{1};
  f32 min_duplex_af{};
  f32 max_duplex_af{};
  f32 min_mapq_mean{};
  std::optional<io::CommandLineInfo> command_line;
};

/**
 * @brief Generate a Panel of Normals VCF from feature TSV files.
 * @param param CLI parameters including file paths, thread count, and filter thresholds.
 * @details Reads feature files in parallel, aggregates variant statistics into a sharded map,
 *          then writes a sites-only VCF via k-way merge across shards.
 */
void GeneratePon(const GeneratePonParam& param);

}  // namespace xoos::svc
