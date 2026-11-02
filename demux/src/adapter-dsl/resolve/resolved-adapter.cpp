#include "adapter-dsl/resolve/resolved-adapter.h"

#include <algorithm>
#include <cstddef>
#include <filesystem>  // NOLINT(build/c++17)
#include <optional>
#include <ranges>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <kseq++/seqio.hpp>

#include <xoos/error/error.h>
#include <xoos/util/sequence-functions.h>

#include "adapter-dsl/model/adapter-model.h"

namespace xoos::demux::dsl {

namespace {

/// Read every sequence from a FASTA file into a list. Throws if the file is unreadable/empty.
std::vector<std::string> ReadFastaPool(const fs::path& path) {
  if (!fs::exists(path)) {
    throw error::Error("Pool FASTA '{}' does not exist", path.string());
  }
  std::vector<std::string> sequences;
  klibpp::SeqStreamIn stream{path.string().c_str()};
  klibpp::KSeq record;
  while (true) {
    stream >> record;
    if (stream.fail()) {
      break;
    }
    sequences.push_back(record.seq);
  }
  if (sequences.empty()) {
    throw error::Error("Pool FASTA '{}' contains no sequences", path.string());
  }
  return sequences;
}

/// Apply reverse complement to every sequence in a pool when the orientation is reverse.
std::vector<std::string> Orient(std::vector<std::string> pool, const Orientation orientation) {
  if (orientation == Orientation::kReverse) {
    for (auto& seq : pool) {
      seq = sequence::ReverseComplement(seq);
    }
  }
  return pool;
}

/// Load an id's forward-strand candidate pool from its FASTA path.
std::vector<std::string> LoadIdPool(const std::string& name,
                                    const IdToken& id,
                                    const AdapterEntry& entry,
                                    const fs::path& resource_root) {
  if (entry.sequence.empty()) {
    throw error::Error("Adapter '{}': id '{}' has no 'sequence' FASTA path", name, id.base);
  }
  return ReadFastaPool(resource_root / entry.sequence);
}

/**
 * Resolve a single leaf token into a runtime element. When @p resource_root is empty the id pool is
 * left unloaded (structure-only mode: see @ref ResolveStructure); otherwise the id's FASTA pool is
 * loaded and oriented.
 */
ResolvedElement ResolveLeaf(const std::string& name,
                            const LeafToken& leaf,
                            const CompiledAdapter& adapter,
                            const std::optional<fs::path>& resource_root) {
  if (const auto* const fixed = std::get_if<FixedToken>(&leaf)) {
    return ResolvedFixed{fixed->sequence};
  }
  if (const auto* const insert = std::get_if<InsertToken>(&leaf)) {
    return ResolvedInsert{insert->orientation};
  }
  const auto& id = std::get<IdToken>(leaf);
  const AdapterEntry* const entry = adapter.Find(id.base);
  if (entry == nullptr) {
    throw error::Error("Adapter '{}': id '{}' has no adapters entry to resolve", name, id.base);
  }
  std::vector<std::string> pool;
  if (resource_root.has_value()) {
    pool = Orient(LoadIdPool(name, id, *entry, *resource_root), id.orientation);
  }
  // bait / read_prefix / read_suffix are derived later by Infer(); initialize them explicitly.
  return ResolvedId{
      id.base, entry->type, id.orientation, std::move(pool), *entry, "", "", "", std::nullopt, std::nullopt};
}

/// Append the resolved elements of one leaf-token sequence in order.
void ResolveInto(const std::string& name,
                 const std::vector<LeafToken>& leaves,
                 const CompiledAdapter& adapter,
                 const std::optional<fs::path>& resource_root,
                 std::vector<ResolvedElement>& out) {
  out.reserve(out.size() + leaves.size());
  for (const LeafToken& leaf : leaves) {
    out.push_back(ResolveLeaf(name, leaf, adapter, resource_root));
  }
}

/**
 * Walk a compiled adapter's read structure into an ordered element list. Pools are loaded when
 * @p resource_root is set and skipped when it is empty (structure-only). This is the single walk
 * shared by @ref ResolveAdapter (pool-loading) and @ref ResolveStructure (structure-only).
 */
ResolvedAdapter WalkStructure(const std::string& name,
                              const CompiledAdapter& adapter,
                              const std::optional<fs::path>& resource_root) {
  ResolvedAdapter resolved;
  resolved.name = name;
  resolved.type = adapter.type;
  resolved.is_duplex = adapter.IsDuplex();
  for (const Token& token : adapter.read_structure) {
    std::visit(
        [&name, &adapter, &resource_root, &resolved]<typename T>(const T& value) {
          if constexpr (std::is_same_v<T, ArmToken>) {
            ResolveInto(name, value.tokens, adapter, resource_root, resolved.elements);
          } else {
            resolved.elements.push_back(ResolveLeaf(name, LeafToken{value}, adapter, resource_root));
          }
        },
        token);
  }
  return resolved;
}

}  // namespace

ResolvedAdapter ResolveAdapter(const std::string& name, const CompiledAdapter& adapter, const fs::path& resource_root) {
  return WalkStructure(name, adapter, resource_root);
}

ResolvedAdapter ResolveStructure(const std::string& name, const CompiledAdapter& adapter) {
  return WalkStructure(name, adapter, std::nullopt);
}

std::vector<ResolvedAdapter> Resolve(const AdapterBundle& bundle, const fs::path& resource_root) {
  std::vector<ResolvedAdapter> resolved;
  resolved.reserve(bundle.adapters.size());
  for (const auto& [name, adapter] : bundle.adapters) {
    resolved.push_back(ResolveAdapter(name, adapter, resource_root));
  }
  return resolved;
}

namespace {

/// A positive length knob as a size, or 0 when unset/non-positive.
std::size_t AsLength(const std::optional<s64>& value) {
  return value.has_value() && value.value() > 0 ? static_cast<std::size_t>(value.value()) : 0;
}

/**
 * The last @p length bases of the ResolvedFixed immediately before @p pos, or empty. When the
 * neighbor is fixed but shorter than @p length the whole fixed is returned (matches v1 baiting).
 */
std::string PrecedingFixedTail(const std::vector<ResolvedElement>& elements,
                               const std::size_t pos,
                               const std::size_t length) {
  if (pos == 0 || length == 0) {
    return {};
  }
  const auto* const fixed = std::get_if<ResolvedFixed>(&elements[pos - 1]);
  if (fixed == nullptr) {
    return {};
  }
  if (fixed->sequence.size() < length) {
    return fixed->sequence;
  }
  return fixed->sequence.substr(fixed->sequence.size() - length);
}

/// The first @p length bases of the ResolvedFixed immediately after @p pos, or empty.
std::string FollowingFixedHead(const std::vector<ResolvedElement>& elements,
                               const std::size_t pos,
                               const std::size_t length) {
  if (pos + 1 >= elements.size() || length == 0) {
    return {};
  }
  const auto* const fixed = std::get_if<ResolvedFixed>(&elements[pos + 1]);
  if (fixed == nullptr) {
    return {};
  }
  return fixed->sequence.substr(0, std::min(length, fixed->sequence.size()));
}

/// Position of the insert nearest to @p pos, or empty when the structure has no insert.
std::optional<std::size_t> NearestInsert(const std::vector<ResolvedElement>& elements, const std::size_t pos) {
  std::optional<std::size_t> nearest;
  std::size_t nearest_distance = elements.size();
  for (const std::size_t index : std::views::iota(std::size_t{0}, elements.size())) {
    if (!std::holds_alternative<ResolvedInsert>(elements[index])) {
      continue;
    }
    const std::size_t distance = index > pos ? index - pos : pos - index;
    if (distance < nearest_distance) {
      nearest = index;
      nearest_distance = distance;
    }
  }
  return nearest;
}

/**
 * Fill @p id's derived-affix fields from the fixeds adjacent to index @p pos. Baits are stored in
 * read orientation (as read); the structural neighbors are stored as element indices (see
 * preceding_element / following_element); UMI prefix/suffix are stored forward-strand (RC-swapped
 * for a reverse id) because the LUT pipeline applies them before the reverse-complement transform.
 */
void InferId(const std::vector<ResolvedElement>& elements, const std::size_t pos, ResolvedId& id) {
  const bool reverse = id.orientation == Orientation::kReverse;

  // Bait side is structural: the fixed outside the SID/insert boundary. Orientation selects only
  // the corresponding length knob; it does not decide which neighboring fixed supplies the bait.
  const std::size_t bait_length = AsLength(reverse ? id.params.bait_length_reverse : id.params.bait_length_forward);
  id.bait.clear();
  const std::optional<std::size_t> insert = NearestInsert(elements, pos);

  // Adjacent-element references, in read-structure order: the index of the neighbor on each side
  // (any element type), or nullopt only at the structure boundary (this id is the first/last element).
  id.preceding_element = pos == 0 ? std::nullopt : std::optional<std::size_t>{pos - 1};
  id.following_element = pos + 1 < elements.size() ? std::optional<std::size_t>{pos + 1} : std::nullopt;

  if (insert.has_value()) {
    id.bait =
        pos < *insert ? PrecedingFixedTail(elements, pos, bait_length) : FollowingFixedHead(elements, pos, bait_length);
  }

  // UMI flank affixes, in read-structure order: prefix = tail of preceding fixed, suffix = head of following.
  const std::string read_prefix = PrecedingFixedTail(elements, pos, AsLength(id.params.additional_prefix_length));
  const std::string read_suffix = FollowingFixedHead(elements, pos, AsLength(id.params.additional_suffix_length));
  if (reverse) {
    id.read_prefix = sequence::ReverseComplement(read_suffix);
    id.read_suffix = sequence::ReverseComplement(read_prefix);
  } else {
    id.read_prefix = read_prefix;
    id.read_suffix = read_suffix;
  }
}

}  // namespace

void Infer(ResolvedAdapter& adapter) {
  for (std::size_t i = 0; i < adapter.elements.size(); ++i) {
    if (auto* const id = std::get_if<ResolvedId>(&adapter.elements[i])) {
      InferId(adapter.elements, i, *id);
    }
  }
}

void Infer(std::vector<ResolvedAdapter>& adapters) {
  for (ResolvedAdapter& adapter : adapters) {
    Infer(adapter);
  }
}

}  // namespace xoos::demux::dsl
