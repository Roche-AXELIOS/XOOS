#pragma once

#include <cstddef>
#include <filesystem>  // NOLINT(build/c++17)
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include <xoos/types/int.h>

#include "adapter-dsl/model/adapter-model.h"

namespace fs = std::filesystem;

namespace xoos::demux::dsl {

/**
 * @brief The runtime view of a compiled adapter: candidate pools loaded and oriented, ready for
 * the matcher.
 *
 * Resolution is the last hop. It takes a validated CompiledAdapter (structure + FASTA paths +
 * params) and:
 *   - loads each id's FASTA pool once,
 *   - materializes each token in read-structure order into a concrete element, applying reverse
 *     complement to the pool for a `-` (reverse) id, and
 *   - carries the v1 params through on the sid/umi elements.
 *
 * The result is an ordered, flat element list (arms flattened, loop preserved by position) that a
 * matcher can walk 5'->3' without re-consulting the DSL or the compiled JSON.
 */

/// @brief A fixed literal element: the sequence is the data itself.
struct ResolvedFixed {
  std::string sequence;
};

/**
 * @brief A resolved barcode/UMI element: an oriented candidate pool plus the v1 search params.
 *
 * The derived-affix fields (@ref bait / @ref read_prefix / @ref read_suffix) are empty after
 * Resolve() and populated by Infer(), which takes them positionally from the neighboring
 * ResolvedFixed element(s) per the id's `bait_length` / `additional_prefix_length` /
 * `additional_suffix_length` knobs. Downstream bundles consume these instead of re-deriving from
 * role tags.
 */
struct ResolvedId {
  /// The id (adapters key).
  std::string base;
  /// sid or umi.
  EntryType type{EntryType::kSid};
  Orientation orientation{Orientation::kForward};
  /// Candidate sequences, reverse-complemented if orientation is reverse.
  std::vector<std::string> pool;
  /// The originating entry (max_edit_distance, wiggle, ...).
  AdapterEntry params;

  /**
   * Derived by Infer(): the bait seed, in read orientation (as read, not reverse-complemented). Its source side is
   * structural: an id before its nearest insert uses the tail of its preceding fixed, while an id
   * after its nearest insert uses the head of its following fixed. The id orientation selects the
   * matching `bait_length+` or `bait_length-` knob only. Empty when the id declares no bait length.
   * (A whole trailing fixed used as a bait -- e.g. YS -- is not stored here; it stays a first-class
   * ResolvedFixed the bundle reads positionally.)
   */
  std::string bait;
  /**
   * Derived by Infer(): the FORWARD-STRAND UMI prefix applied to each pool candidate before the
   * reverse-complement transform. Forward id: last `additional_prefix_length` bases of the fixed 5'
   * of the id. Reverse id: reverse complement of the first `additional_suffix_length` bases of the
   * fixed 3' of the id (read-order prefix/suffix swap under RC). Empty when unused.
   */
  std::string read_prefix;
  /**
   * Derived by Infer(): the FORWARD-STRAND UMI suffix applied to each pool candidate before the
   * reverse-complement transform. Forward id: first `additional_suffix_length` bases of the fixed
   * 3' of the id. Reverse id: reverse complement of the last `additional_prefix_length` bases of
   * the fixed 5' of the id. Empty when unused.
   */
  std::string read_suffix;
  /**
   * Index into @ref ResolvedAdapter::elements of the element immediately before this id in read
   * order, or nullopt when this id is the first element. A general structural back-reference: the
   * neighbor may be any element type, and consumers interpret it (e.g. the assembler reads it as a
   * flank only when it is a ResolvedFixed). Referencing the neighbor by position, rather than
   * copying its data, keeps the id anchored to the read structure.
   */
  std::optional<std::size_t> preceding_element;
  /**
   * Index into @ref ResolvedAdapter::elements of the element immediately after this id in read
   * order, or nullopt when this id is the last element. See @ref preceding_element.
   */
  std::optional<std::size_t> following_element;
};

/// @brief The insert element: the payload region, self-describing (no pool).
struct ResolvedInsert {
  std::optional<Orientation> orientation;
};

/// @brief One resolved element in read-structure order.
using ResolvedElement = std::variant<ResolvedFixed, ResolvedId, ResolvedInsert>;

/// @brief A fully resolved adapter: type plus the ordered runtime element list.
struct ResolvedAdapter {
  std::string name;
  DesignType type{DesignType::kSimplex};
  bool is_duplex{false};
  std::vector<ResolvedElement> elements;
};

/**
 * @brief Resolve every adapter in @p bundle against FASTA pools rooted at @p resource_root.
 *
 * Each entry's `sequence` field is a path relative to @p resource_root. Pools are loaded once per
 * (id) and reverse-complemented when the referencing token is reverse.
 *
 * The sample-sheet SID pool is substituted downstream at LUT-generation time (unchanged); it is
 * not an input here. Resolution is purely config/FASTA-driven.
 *
 * @throws xoos::Error if a FASTA is missing/empty or an id has no sequence path.
 */
std::vector<ResolvedAdapter> Resolve(  // NOSONAR(S5536) -- exercised by tests (not in sonar.sources)
    const AdapterBundle& bundle,
    const fs::path& resource_root);

/// @brief Resolve a single named adapter from @p bundle.
ResolvedAdapter ResolveAdapter(const std::string& name, const CompiledAdapter& adapter, const fs::path& resource_root);

/**
 * @brief Build the structural view of a single adapter WITHOUT loading any FASTA pools.
 *
 * Produces the same ordered element list as @ref ResolveAdapter (arms flattened, loop preserved by
 * position, orientation carried) but leaves every @ref ResolvedId::pool empty. This is the input the
 * design assembler needs: derivation (@ref Infer) reads only the neighboring ResolvedFixed literals
 * and the id knobs, never the pools, and the assembler emits FASTA *paths* into the blueprint rather
 * than sequences. Because it touches no filesystem it runs on the VFS runtime path.
 */
ResolvedAdapter ResolveStructure(const std::string& name, const CompiledAdapter& adapter);

/**
 * @brief Annotate every ResolvedId in @p adapter with its derived affixes, in place.
 *
 * Infer() walks the ordered element list and, for each id, reads the neighboring ResolvedFixed
 * element(s) to fill @ref ResolvedId::bait / @ref ResolvedId::read_prefix / @ref
 * ResolvedId::read_suffix from the id's `bait_length` / `additional_prefix_length` /
 * `additional_suffix_length` knobs. This is the only step that derives data from the config; it is
 * positional (no role tags) and idempotent. Ids that declare no length knobs are left untouched.
 *
 * Runs after Resolve() and before resource construction.
 */
void Infer(ResolvedAdapter& adapter);

/// @brief Run @ref Infer over every adapter in @p adapters, in place.
void Infer(std::vector<ResolvedAdapter>& adapters);  // NOSONAR(S5536) -- exercised by tests (not in sonar.sources)

}  // namespace xoos::demux::dsl
