#include "phase-block-builder.h"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <optional>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <xoos/log/logging.h>
#include <xoos/types/int.h>

namespace xoos::svc::phaser {
namespace {

// A single phased pair is the smallest useful phase block.
constexpr size_t kMinVariantsToPhase = 2;

// Disjoint-set structure used to collapse pairwise phase edges into connected components.
class UnionFind {
 public:
  // Start with one component per variant index.
  explicit UnionFind(const s32 count) : _parent(static_cast<size_t>(count)), _rank(static_cast<size_t>(count), 0) {
    std::ranges::iota(_parent, 0);
  }

  // Find a component representative and compress the traversed path.
  s32 Find(s32 value) {
    s32 root = value;
    while (_parent[static_cast<size_t>(root)] != root) {
      root = _parent[static_cast<size_t>(root)];
    }
    while (_parent[static_cast<size_t>(value)] != root) {
      const s32 next = _parent[static_cast<size_t>(value)];
      _parent[static_cast<size_t>(value)] = root;
      value = next;
    }
    return root;
  }

  // Merge two components, keeping the shallower tree as a child.
  void Unite(const s32 lhs, const s32 rhs) {
    s32 lhs_root = Find(lhs);
    s32 rhs_root = Find(rhs);
    if (lhs_root == rhs_root) {
      return;
    }
    if (_rank[static_cast<size_t>(lhs_root)] < _rank[static_cast<size_t>(rhs_root)]) {
      std::swap(lhs_root, rhs_root);
    }
    _parent[static_cast<size_t>(rhs_root)] = lhs_root;
    if (_rank[static_cast<size_t>(lhs_root)] == _rank[static_cast<size_t>(rhs_root)]) {
      ++_rank[static_cast<size_t>(lhs_root)];
    }
  }

 private:
  // Parent links and ranks are indexed by the original variant index.
  vec<s32> _parent;
  vec<s32> _rank;
};

// For each variant, store neighboring variants and whether their edge is cis.
using PhaseAdjacency = std::unordered_map<s32, vec<std::pair<s32, bool>>>;

std::unordered_map<s32, vec<s32>> BuildPhaseComponents(const vec<PhasingResult>& results, const size_t variant_count) {
  UnionFind union_find(static_cast<s32>(variant_count));
  std::unordered_set<s32> phased_variants;

  // Join every phased pair; unphased results do not contribute to a block.
  for (const auto& result : results) {
    if (result.is_phased) {
      union_find.Unite(result.variant_idx_a, result.variant_idx_b);
      phased_variants.insert(result.variant_idx_a);
      phased_variants.insert(result.variant_idx_b);
    }
  }

  // Keep only variants that occur in a phased edge and group them by representative.
  std::unordered_map<s32, vec<s32>> components;
  for (const s32 variant_index : phased_variants) {
    components[union_find.Find(variant_index)].push_back(variant_index);
  }
  return components;
}

PhaseAdjacency BuildPhaseAdjacency(const vec<PhasingResult>& results) {
  PhaseAdjacency adjacency;

  // Store edges in both directions so the coloring traversal can start anywhere.
  for (const auto& result : results) {
    if (result.is_phased) {
      adjacency[result.variant_idx_a].emplace_back(result.variant_idx_b, result.is_cis);
      adjacency[result.variant_idx_b].emplace_back(result.variant_idx_a, result.is_cis);
    }
  }
  return adjacency;
}

std::optional<std::unordered_map<s32, s32>> TwoColorComponent(const vec<s32>& indices,
                                                              const PhaseAdjacency& adjacency) {
  std::unordered_map<s32, s32> color;
  for (const s32 variant_index : indices) {
    color[variant_index] = -1;
  }

  // A cis edge keeps a haplotype color; a trans edge flips it.
  std::queue<s32> bfs_queue;
  color[indices[0]] = 0;
  bfs_queue.push(indices[0]);

  while (!bfs_queue.empty()) {
    const s32 current = bfs_queue.front();
    bfs_queue.pop();

    const auto adjacent = adjacency.find(current);
    if (adjacent == adjacency.end()) {
      continue;
    }

    for (const auto& [neighbor, is_cis] : adjacent->second) {
      if (!color.contains(neighbor)) {
        continue;
      }
      const s32 expected_color = is_cis ? color[current] : (1 - color[current]);
      if (color[neighbor] == -1) {
        color[neighbor] = expected_color;
        bfs_queue.push(neighbor);
      } else if (color[neighbor] != expected_color) {
        // Conflicting evidence means no consistent two-haplotype assignment exists.
        return std::nullopt;
      }
    }
  }
  return color;
}

}  // namespace

vec<PhaseBlock> BuildPhaseBlocks(const vec<PhasingResult>& results, const vec<Variant>& variants) {
  // Build both the connectivity view and the edge view used for haplotype coloring.
  auto components = BuildPhaseComponents(results, variants.size());
  const PhaseAdjacency adjacency = BuildPhaseAdjacency(results);

  vec<PhaseBlock> blocks;
  s32 coloring_conflicts = 0;
  for (auto& [root, indices] : components) {
    // Singleton components cannot represent a phased relationship.
    if (indices.size() < kMinVariantsToPhase) {
      continue;
    }

    // Make output deterministic and remove any duplicate index introduced by the graph traversal.
    std::ranges::sort(indices, [&variants](const s32 lhs, const s32 rhs) {
      return variants[static_cast<size_t>(lhs)].pos < variants[static_cast<size_t>(rhs)].pos;
    });
    const auto duplicates = std::ranges::unique(indices);
    indices.erase(duplicates.begin(), duplicates.end());

    const auto color = TwoColorComponent(indices, adjacency);
    PhaseBlock block;
    block.variant_indices = indices;
    if (color.has_value()) {
      // Store one color per sorted variant for downstream PGT/PS annotation.
      block.haplotypes.resize(indices.size());
      for (size_t i = 0; i < indices.size(); ++i) {
        block.haplotypes[i] = color->at(indices[i]);
      }
    } else {
      // Preserve the block for phase-set output, but omit unreliable haplotypes.
      ++coloring_conflicts;
    }
    blocks.push_back(std::move(block));
  }

  if (coloring_conflicts > 0) {
    Logging::Warn("{} phase blocks had conflicting cis/trans evidence — haplotype assignment skipped",
                  coloring_conflicts);
  }
  return blocks;
}

}  // namespace xoos::svc::phaser
