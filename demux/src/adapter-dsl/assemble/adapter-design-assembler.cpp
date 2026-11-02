#include "adapter-dsl/assemble/adapter-design-assembler.h"

#include <cstddef>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>
#include <xoos/util/sequence-functions.h>

#include "adapter-design/adapter-component-indices.h"
#include "adapter-dsl/resolve/resolved-adapter.h"

namespace xoos::demux::dsl {

namespace {

// --- small field helpers -------------------------------------------------------------------------

SeqTransform TransformOf(const Orientation orientation) {
  return orientation == Orientation::kReverse ? SeqTransform::kReverseComplement : SeqTransform::kNone;
}

u32 EntryEdist(const AdapterEntry& entry) {
  return static_cast<u32>(entry.max_edit_distance.value_or(0));
}

s32 WiggleLeft(const std::optional<Wiggle>& w) {
  return w.has_value() ? static_cast<s32>(w->left) : 0;
}

s32 WiggleRight(const std::optional<Wiggle>& w) {
  return w.has_value() ? static_cast<s32>(w->right) : 0;
}

const std::optional<Wiggle>& WiggleFor(const AdapterEntry& entry, const Orientation orientation) {
  return orientation == Orientation::kReverse ? entry.wiggle_reverse : entry.wiggle_forward;
}

SearchWindow MakeSearch(const s32 left, const s32 right) {
  return SearchWindow{left, right};
}

const AdapterEntry* FindEntry(const CompiledAdapter& adapter, const std::string& base) {
  return adapter.Find(base);
}

const AdapterEntry& RequireEntry(const std::string& name, const CompiledAdapter& adapter, const std::string& base) {
  const AdapterEntry* const entry = FindEntry(adapter, base);
  if (entry == nullptr) {
    throw error::Error("Adapter '{}': missing required id '{}'", name, base);
  }
  return *entry;
}

const AdapterEntry& RequireEntryOfType(const std::string& name,
                                       const CompiledAdapter& adapter,
                                       const std::string& base,
                                       const EntryType type) {
  const AdapterEntry& entry = RequireEntry(name, adapter, base);
  if (entry.type != type) {
    throw error::Error("Adapter '{}': id '{}' has incompatible adapter entry type", name, base);
  }
  return entry;
}

// --- component factories -------------------------------------------------------------------------

/**
 * A fixed (single inline literal) component. The trimmer distinguishes fixed literals by position,
 * not by a sub-role tag.
 */
FixedComponentDefinition MakeFixed(const std::string& sequence) {
  return FixedComponentDefinition{sequence};
}

/**
 * A barcode-pool (multi-record FASTA) component for a sid/umi id at a given orientation. Optional
 * @p prefix / @p suffix are forward-strand literal affixes (applied before the transform by the LUT
 * pipeline); they carry UMI flank extensions and are supplied already forward-strand by Infer.
 */
MatcherComponentDefinition MakePool(const MatcherKind kind,
                                    const AdapterEntry& entry,
                                    const Orientation orientation,
                                    const std::string_view prefix,
                                    const std::string_view suffix) {
  MatcherComponentDefinition bd;
  bd.kind = kind;
  bd.sequences = entry.sequence;
  bd.max_edit_distance = EntryEdist(entry);
  bd.transform = TransformOf(orientation);
  bd.prefix = prefix;
  bd.suffix = suffix;
  const std::optional<Wiggle>& w = WiggleFor(entry, orientation);
  bd.search = MakeSearch(WiggleLeft(w), WiggleRight(w));
  return bd;
}

/// A sid pool component with no affixes.
MatcherComponentDefinition MakeSid(const AdapterEntry& entry, const Orientation orientation) {
  return MakePool(MatcherKind::kSid, entry, orientation, {}, {});
}

/**
 * A SID pool component extended with its adjacent read-order spacer. Matcher affixes are applied
 * before the sequence transform, so reverse-oriented components must reverse-complement and swap
 * the read-order affixes first. This keeps spacer placement structural: the first YS SID consumes
 * its following fixed and the last SID consumes its preceding fixed, regardless of orientation.
 */
MatcherComponentDefinition MakeSidWithReadAffixes(const AdapterEntry& entry,
                                                  const Orientation orientation,
                                                  const std::string& read_prefix,
                                                  const std::string& read_suffix) {
  if (orientation == Orientation::kReverse) {
    return MakePool(MatcherKind::kSid,
                    entry,
                    orientation,
                    sequence::ReverseComplement(read_suffix),
                    sequence::ReverseComplement(read_prefix));
  }
  return MakePool(MatcherKind::kSid, entry, orientation, read_prefix, read_suffix);
}

// --- derivation lookup (single derivation site is ResolveStructure + Infer) ----------------------

/**
 * The Infer'd ResolvedId for id @p base at @p orientation, or throws. This is how the assembler
 * consumes derived values (UMI flank affixes, sid baits) without re-deriving them: Infer computed
 * them once from the resolved structure, and the assembler only reads them back.
 */
const ResolvedId& FindResolvedId(const std::string& name,
                                 const ResolvedAdapter& resolved,
                                 const std::string& base,
                                 const Orientation orientation) {
  for (const ResolvedElement& element : resolved.elements) {
    if (const auto* const id = std::get_if<ResolvedId>(&element)) {
      if (id->base == base && id->orientation == orientation) {
        return *id;
      }
    }
  }
  throw error::Error("Adapter '{}': no resolved id '{}' in the required orientation", name, base);
}

/// All resolved ids named @p base, in read-structure order. Throws unless exactly @p count are present.
vec<const ResolvedId*> FindResolvedIds(const std::string& name,
                                       const ResolvedAdapter& resolved,
                                       const std::string& base,
                                       const std::size_t count) {
  vec<const ResolvedId*> ids;
  for (const ResolvedElement& element : resolved.elements) {
    if (const auto* const id = std::get_if<ResolvedId>(&element); id != nullptr && id->base == base) {
      ids.push_back(id);
    }
  }
  if (ids.size() != count) {
    throw error::Error("Adapter '{}': expected {} resolved '{}' ids, found {}", name, count, base, ids.size());
  }
  return ids;
}

/// A UMI pool component whose flank affixes are the Infer'd (already forward-strand) values.
MatcherComponentDefinition MakeUmi(const ResolvedId& umi) {
  return MakePool(MatcherKind::kUmi, umi.params, umi.orientation, umi.read_prefix, umi.read_suffix);
}

// --- read-structure walkers (geometry only; no value derivation) ---------------------------------

const FixedToken* AsFixed(const LeafToken& leaf) {
  return std::get_if<FixedToken>(&leaf);
}

const IdToken* AsId(const LeafToken& leaf) {
  return std::get_if<IdToken>(&leaf);
}

/**
 * Flatten the read structure to its leaf tokens in read order (arms expanded in place). Used by the
 * flat-strand families (YSU / YS / simplex / simplex-10x), which have no arms.
 */
vec<LeafToken> FlattenLeaves(const vec<Token>& read_structure) {
  vec<LeafToken> leaves;
  for (const Token& token : read_structure) {
    std::visit(
        [&leaves]<typename T>(const T& value) {
          if constexpr (std::is_same_v<T, ArmToken>) {
            for (const LeafToken& leaf : value.tokens) {
              leaves.emplace_back(leaf);
            }
          } else {
            leaves.emplace_back(value);
          }
        },
        token);
  }
  return leaves;
}

/// Positions of all id tokens named @p base, in read-structure order. Throws unless @p count are present.
vec<std::size_t> FindIdIndices(const std::string& name,
                               const vec<LeafToken>& leaves,
                               const std::string& base,
                               const std::size_t count) {
  vec<std::size_t> indices;
  for (std::size_t i = 0; i < leaves.size(); ++i) {
    if (const IdToken* const id = AsId(leaves[i])) {
      if (id->base == base) {
        indices.push_back(i);
      }
    }
  }
  if (indices.size() != count) {
    throw error::Error("Adapter '{}': expected {} '{}' id tokens, found {}", name, count, base, indices.size());
  }
  return indices;
}

/// Orientation of the id token at @p index. The caller obtained the index from FindIdIndices.
Orientation IdOrientationAt(const vec<LeafToken>& leaves, const std::size_t index) {
  return AsId(leaves[index])->orientation;
}

/// Position of the first id token with @p base within @p leaves, or throws.
std::size_t FindIdIndex(const std::string& name, const vec<LeafToken>& leaves, const std::string& base) {
  const vec<std::size_t> indices = FindIdIndices(name, leaves, base, 1);
  return indices.front();
}

/// The first two positions of id token @p base, in read-structure order.
std::pair<std::size_t, std::size_t> FindIdPair(const std::string& name,
                                               const vec<LeafToken>& leaves,
                                               const std::string& base) {
  const vec<std::size_t> indices = FindIdIndices(name, leaves, base, 2);
  return {indices.front(), indices.back()};
}

/// The first two resolved ids named @p base, in read-structure order.
std::pair<const ResolvedId&, const ResolvedId&> FindResolvedIdPair(const std::string& name,
                                                                   const ResolvedAdapter& resolved,
                                                                   const std::string& base) {
  const vec<const ResolvedId*> ids = FindResolvedIds(name, resolved, base, 2);
  return {*ids.front(), *ids.back()};
}

/// @brief Look up an id token at an index previously returned by FindIdIndices.
const IdToken& IdAt(const vec<LeafToken>& leaves, const std::size_t index) {
  const IdToken* const id = AsId(leaves[index]);
  if (id == nullptr) {
    throw error::Error("Internal adapter assembly error: expected an id at position {}", index);
  }
  return *id;
}

std::size_t FindInsertIndex(const std::string& name, const vec<LeafToken>& leaves) {
  for (std::size_t i = 0; i < leaves.size(); ++i) {
    if (std::holds_alternative<InsertToken>(leaves[i])) {
      return i;
    }
  }
  throw error::Error("Adapter '{}': no insert token", name);
}

bool IsSidId(const CompiledAdapter& adapter, const IdToken& id) {
  const AdapterEntry* const entry = FindEntry(adapter, id.base);
  return entry != nullptr && entry->type == EntryType::kSid;
}

std::size_t FindSidIndexOnInsertSide(const std::string& name,
                                     const CompiledAdapter& adapter,
                                     const vec<LeafToken>& leaves,
                                     const std::size_t insert_idx,
                                     const bool before_insert) {
  std::optional<std::size_t> result;
  for (const auto& [i, leaf] : std::views::enumerate(leaves)) {
    const auto index = static_cast<std::size_t>(i);
    if (before_insert != (index < insert_idx)) {
      continue;
    }
    const IdToken* const id = AsId(leaf);
    if (id != nullptr && IsSidId(adapter, *id)) {
      if (result.has_value()) {
        throw error::Error("Adapter '{}': simplex structure has multiple SID ids on one side of the insert", name);
      }
      result = index;
    }
  }
  if (!result.has_value()) {
    throw error::Error("Adapter '{}': simplex structure is missing a SID id on one side of the insert", name);
  }
  return *result;
}

/// The full fixed literal at leaf index @p idx, or throws if the leaf is missing or not fixed.
std::string RequireFixedAt(const std::string& name, const vec<LeafToken>& leaves, const std::size_t idx) {
  if (idx >= leaves.size()) {
    throw error::Error("Adapter '{}': expected a fixed literal at position {}", name, idx);
  }
  const FixedToken* const fixed = AsFixed(leaves[idx]);
  if (fixed == nullptr) {
    throw error::Error("Adapter '{}': expected a fixed literal at position {}", name, idx);
  }
  return fixed->sequence;
}

// --- duplex arm walkers --------------------------------------------------------------------------

const ArmToken& RequireArm(const std::string& name, const CompiledAdapter& adapter, const std::size_t which) {
  std::size_t seen = 0;
  for (const Token& token : adapter.read_structure) {
    if (const auto* const arm = std::get_if<ArmToken>(&token)) {
      if (seen == which) {
        return *arm;
      }
      ++seen;
    }
  }
  throw error::Error("Adapter '{}': missing duplex arm #{}", name, which + 1);
}

/// Leading top-level fixed (before the first arm) = runway.
std::string DuplexRunway(const std::string& name, const CompiledAdapter& adapter) {
  for (const Token& token : adapter.read_structure) {
    if (const auto* const fixed = std::get_if<FixedToken>(&token)) {
      return fixed->sequence;
    }
    if (std::holds_alternative<ArmToken>(token)) {
      break;
    }
  }
  throw error::Error("Adapter '{}': no leading runway fixed", name);
}

/// The single top-level fixed between the two arms = loop (identified by position).
std::string DuplexLoop(const std::string& name, const CompiledAdapter& adapter) {
  bool after_first_arm = false;
  for (const Token& token : adapter.read_structure) {
    if (std::holds_alternative<ArmToken>(token)) {
      if (after_first_arm) {
        break;
      }
      after_first_arm = true;
      continue;
    }
    if (after_first_arm) {
      if (const auto* const fixed = std::get_if<FixedToken>(&token)) {
        return fixed->sequence;
      }
    }
  }
  throw error::Error("Adapter '{}': no loop fixed between its arms", name);
}

/// Outermost fixed in an arm: first fixed when @p from_front, else last fixed.
std::string ArmOuterFixed(const std::string& name, const ArmToken& arm, const bool from_front) {
  const FixedToken* found = nullptr;
  for (const LeafToken& leaf : arm.tokens) {
    if (const FixedToken* const fixed = AsFixed(leaf)) {
      found = fixed;
      if (from_front) {
        break;
      }
    }
  }
  if (found == nullptr) {
    throw error::Error("Adapter '{}': duplex arm has no anchor fixed", name);
  }
  return found->sequence;
}

/// Number of fixed leaves in an arm.
std::size_t ArmFixedCount(const ArmToken& arm) {
  std::size_t count = 0;
  for (const LeafToken& leaf : arm.tokens) {
    if (AsFixed(leaf) != nullptr) {
      ++count;
    }
  }
  return count;
}

/// Orientation of the first id token with @p base within an arm.
Orientation ArmIdOrientation(const std::string& name, const ArmToken& arm, const std::string& base) {
  return IdAt(arm.tokens, FindIdIndex(name, arm.tokens, base)).orientation;
}

// --- per-family assembly -------------------------------------------------------------------------
//
// Geometry (component read order, duplex arm layout) comes from read_structure. Derived
// values (UMI affixes, YSU sid baits) come from the Infer'd ResolvedAdapter -- never re-derived here.

// SBX-D (duplex). Read order: runway, sid+, anchor, loop, anchor, sid-.
AdapterDesign BuildDuplex(const std::string& name, const CompiledAdapter& adapter) {
  const ArmToken& arm1 = RequireArm(name, adapter, 0);
  const ArmToken& arm2 = RequireArm(name, adapter, 1);
  const std::string runway = DuplexRunway(name, adapter);
  const std::string loop = DuplexLoop(name, adapter);
  const std::string anchor_5p = ArmOuterFixed(name, arm1, /*from_front=*/true);
  const std::string anchor_3p = ArmOuterFixed(name, arm2, /*from_front=*/false);
  const AdapterEntry& sid = RequireEntry(name, adapter, "sid");
  const Orientation sid_5p_orient = ArmIdOrientation(name, arm1, "sid");
  const Orientation sid_3p_orient = ArmIdOrientation(name, arm2, "sid");

  AdapterDesign design;
  design.name = name;
  design.type = AdapterType::kDuplex;

  // Single components vector in molecule read order: runway, sid+, anchor, loop, anchor, sid-.
  design.components.emplace_back(MakeFixed(runway));
  design.components.emplace_back(MakeSid(sid, sid_5p_orient));
  design.components.emplace_back(MakeFixed(anchor_5p));
  design.components.emplace_back(MakeFixed(loop));
  design.components.emplace_back(MakeFixed(anchor_3p));
  design.components.emplace_back(MakeSid(sid, sid_3p_orient));
  return design;
}

// duplex_stem. Read order: runway, sid+, anchor, stem, loop, sid-, anchor, stem.
AdapterDesign BuildDuplexStem(const std::string& name, const CompiledAdapter& adapter) {
  const ArmToken& arm1 = RequireArm(name, adapter, 0);
  const ArmToken& arm2 = RequireArm(name, adapter, 1);
  // Each arm must carry both an anchor and a stem fixed component.
  constexpr std::size_t kMinArmFixedCount = 2;
  if (ArmFixedCount(arm1) < kMinArmFixedCount || ArmFixedCount(arm2) < kMinArmFixedCount) {
    throw error::Error("Adapter '{}': duplex_stem requires an anchor and a stem fixed in each arm", name);
  }
  const std::string runway = DuplexRunway(name, adapter);
  const std::string loop = DuplexLoop(name, adapter);
  const std::string anchor_5p = ArmOuterFixed(name, arm1, /*from_front=*/true);
  const std::string stem_5p = ArmOuterFixed(name, arm1, /*from_front=*/false);
  const std::string stem_3p = ArmOuterFixed(name, arm2, /*from_front=*/true);
  const std::string anchor_3p = ArmOuterFixed(name, arm2, /*from_front=*/false);
  const AdapterEntry& sid = RequireEntry(name, adapter, "sid");
  const Orientation sid_5p_orient = ArmIdOrientation(name, arm1, "sid");
  const Orientation sid_3p_orient = ArmIdOrientation(name, arm2, "sid");

  AdapterDesign design;
  design.name = name;
  design.type = AdapterType::kDuplexStem;

  // Component order: runway, sid+, anchor, stem, loop, sid-, anchor, stem.
  design.components.emplace_back(MakeFixed(runway));
  design.components.emplace_back(MakeSid(sid, sid_5p_orient));
  design.components.emplace_back(MakeFixed(anchor_5p));
  design.components.emplace_back(MakeFixed(stem_5p));
  design.components.emplace_back(MakeFixed(loop));
  design.components.emplace_back(MakeSid(sid, sid_3p_orient));
  design.components.emplace_back(MakeFixed(anchor_3p));
  design.components.emplace_back(MakeFixed(stem_3p));
  return design;
}

// duplex_umi. duplex_stem components plus {umi 5', umi 3'} (umi1, adjacent to the anchors). UMI
// flank affixes come from the Infer'd resolved ids.
AdapterDesign BuildDuplexUmi(const std::string& name, const CompiledAdapter& adapter, const ResolvedAdapter& resolved) {
  // Start from the stem layout [runway, sid+, anchor, stem, loop, sid-, anchor, stem] and insert the
  // two UMIs at the tail of each arm (after the 5' stem, and after the 3' stem).
  AdapterDesign design = BuildDuplexStem(name, adapter);
  design.type = AdapterType::kDuplexUMI;

  const ResolvedId& umi_5p = FindResolvedId(name, resolved, "umi1", Orientation::kForward);
  const ResolvedId& umi_3p = FindResolvedId(name, resolved, "umi1", Orientation::kReverse);

  // Read order becomes: runway, sid+, anchor, stem, loop, umi1+, sid-, anchor, stem, umi1-.
  // The 5' UMI is inserted before the reverse SID; the 3' UMI goes at the very end.
  design.components.emplace(design.components.begin() + adapter_component_index::DuplexStem::kSidRev, MakeUmi(umi_5p));
  design.components.emplace_back(MakeUmi(umi_3p));
  return design;
}

AdapterDesign BuildYsu(const std::string& name, const CompiledAdapter& adapter, const ResolvedAdapter& resolved) {
  using enum Orientation;
  const AdapterEntry& sid = RequireEntry(name, adapter, "sid");

  // umi1/umi2 presence is validated by the FindResolvedId lookups below.
  const auto [sid_first, sid_last] = FindResolvedIdPair(name, resolved, "sid");
  const ResolvedId& umi_5p = FindResolvedId(name, resolved, "umi1", kForward);
  const ResolvedId& umi_3p = FindResolvedId(name, resolved, "umi2", kReverse);

  AdapterDesign design;
  design.name = name;
  design.type = AdapterType::kYsu;

  MatcherComponentDefinition first_sid = MakeSid(sid, sid_first.orientation);
  first_sid.bait = sid_first.bait;
  MatcherComponentDefinition last_sid = MakeSid(sid, sid_last.orientation);
  last_sid.bait = sid_last.bait;

  design.components.emplace_back(std::move(first_sid));
  design.components.emplace_back(MakeUmi(umi_5p));
  design.components.emplace_back(MakeUmi(umi_3p));
  design.components.emplace_back(std::move(last_sid));
  return design;
}

AdapterDesign BuildYs(const std::string& name, const CompiledAdapter& adapter) {
  const vec<LeafToken> leaves = FlattenLeaves(adapter.read_structure);
  const AdapterEntry& sid = RequireEntry(name, adapter, "sid");
  const auto [sid_first, sid_last] = FindIdPair(name, leaves, "sid");

  const std::string runway_5p = RequireFixedAt(name, leaves, 0);
  const std::string spacer_5p = RequireFixedAt(name, leaves, sid_first + 1);
  const std::string spacer_3p = RequireFixedAt(name, leaves, sid_last - 1);
  const std::string bait_3p = RequireFixedAt(name, leaves, sid_last + 1);

  AdapterDesign design;
  design.name = name;
  design.type = AdapterType::kYs;

  design.components.emplace_back(MakeFixed(runway_5p));
  design.components.emplace_back(MakeSidWithReadAffixes(sid, IdOrientationAt(leaves, sid_first), {}, spacer_5p));
  design.components.emplace_back(MakeFixed(spacer_5p));
  design.components.emplace_back(MakeSidWithReadAffixes(sid, IdOrientationAt(leaves, sid_last), spacer_3p, {}));
  design.components.emplace_back(MakeFixed(bait_3p));
  return design;
}

AdapterDesign BuildSimplex(const std::string& name, const CompiledAdapter& adapter) {
  const vec<LeafToken> leaves = FlattenLeaves(adapter.read_structure);
  const std::size_t insert_idx = FindInsertIndex(name, leaves);
  const std::size_t sid_first = FindSidIndexOnInsertSide(name, adapter, leaves, insert_idx, /*before_insert=*/true);
  const std::size_t sid_last = FindSidIndexOnInsertSide(name, adapter, leaves, insert_idx, /*before_insert=*/false);
  const IdToken& sid_5p_token = *AsId(leaves.at(sid_first));
  const IdToken& sid_3p_token = *AsId(leaves.at(sid_last));
  const AdapterEntry& sid_5p_entry = RequireEntryOfType(name, adapter, sid_5p_token.base, EntryType::kSid);
  const AdapterEntry& sid_3p_entry = RequireEntryOfType(name, adapter, sid_3p_token.base, EntryType::kSid);

  const std::string anchor_5p = RequireFixedAt(name, leaves, 0);
  const std::string stem_5p = RequireFixedAt(name, leaves, sid_first + 1);
  const std::string stem_3p = RequireFixedAt(name, leaves, sid_last - 1);
  const std::string anchor_3p = RequireFixedAt(name, leaves, sid_last + 1);

  AdapterDesign design;
  design.name = name;
  design.type = AdapterType::kSimplex;

  design.components.emplace_back(MakeFixed(anchor_5p));
  design.components.emplace_back(MakeSid(sid_5p_entry, sid_5p_token.orientation));
  design.components.emplace_back(MakeFixed(stem_5p));
  design.components.emplace_back(MakeFixed(stem_3p));
  design.components.emplace_back(MakeSid(sid_3p_entry, sid_3p_token.orientation));
  design.components.emplace_back(MakeFixed(anchor_3p));
  return design;
}

// simplex-10x
AdapterDesign BuildSimplex10x(const std::string& name, const CompiledAdapter& adapter) {
  const vec<LeafToken> leaves = FlattenLeaves(adapter.read_structure);
  const AdapterEntry& sid = RequireEntry(name, adapter, "sid");
  const std::size_t sid_index = FindIdIndex(name, leaves, "sid");

  const std::string anchor_5p = RequireFixedAt(name, leaves, 0);
  const std::string stem_5p = RequireFixedAt(name, leaves, sid_index + 1);

  AdapterDesign design;
  design.name = name;
  design.type = AdapterType::kSimplex10x;

  design.components.emplace_back(MakeFixed(anchor_5p));
  design.components.emplace_back(MakeSid(sid, IdOrientationAt(leaves, sid_index)));
  design.components.emplace_back(MakeFixed(stem_5p));
  return design;
}

// Dispatch to the per-type builder for an already-resolved adapter.
AdapterDesign BuildForType(const std::string& name, const CompiledAdapter& adapter, const ResolvedAdapter& resolved) {
  switch (adapter.type) {
    using enum DesignType;
    case kDuplex:
      return BuildDuplex(name, adapter);
    case kDuplexStem:
      return BuildDuplexStem(name, adapter);
    case kDuplexUmi:
      return BuildDuplexUmi(name, adapter, resolved);
    case kYsu:
      return BuildYsu(name, adapter, resolved);
    case kYs:
      return BuildYs(name, adapter);
    case kSimplex:
      return BuildSimplex(name, adapter);
    case kSimplex10x:
      return BuildSimplex10x(name, adapter);
    default:
      throw error::Error("Unhandled adapter design type '{}' for '{}'", ToString(adapter.type), name);
  }
}

/// Sequence of the ResolvedFixed at @p element, or empty when the id has no fixed neighbor there.
std::string StructuralFixedAt(const vec<ResolvedElement>& elements, const std::optional<std::size_t>& element) {
  if (!element.has_value()) {
    return {};
  }
  const auto* const fixed = std::get_if<ResolvedFixed>(&elements[*element]);
  return fixed == nullptr ? std::string{} : fixed->sequence;
}

/**
 * @brief Carry DSL-derived structural SID context onto the assembled runtime components without
 * changing which fixed literals are consumed by each adapter-specific trimmer.
 *
 * @param[in,out] design   Assembled adapter whose SID matcher components receive structural flanks.
 * @param resolved         Resolved DSL adapter providing the structural read prefix/suffix per SID.
 * @throws error::Error If the number of assembled SID components does not match the resolved SIDs.
 */
void AttachStructuralSidFlanks(AdapterDesign& design, const ResolvedAdapter& resolved) {
  vec<const ResolvedId*> resolved_sids;
  for (const ResolvedElement& element : resolved.elements) {
    if (const auto* const id = std::get_if<ResolvedId>(&element); id != nullptr && id->type == EntryType::kSid) {
      resolved_sids.push_back(id);
    }
  }

  std::size_t sid_index = 0;
  for (AdapterComponentDefinition& component : design.components) {
    auto* const matcher = std::get_if<MatcherComponentDefinition>(&component);
    if (matcher == nullptr || matcher->kind != MatcherKind::kSid) {
      continue;
    }
    if (sid_index >= resolved_sids.size()) {
      throw error::Error("Adapter '{}' assembled more SID components than its resolved DSL", design.name);
    }
    const ResolvedId& id = *resolved_sids[sid_index];
    ++sid_index;
    // TODO: See MatcherComponentDefinition::structural_prefix. Once AdapterDesign carries its
    // ordered elements, the matcher can hold a neighbor index instead of these materialized copies.
    matcher->structural_prefix = StructuralFixedAt(resolved.elements, id.preceding_element);
    matcher->structural_suffix = StructuralFixedAt(resolved.elements, id.following_element);
  }
  if (sid_index != resolved_sids.size()) {
    throw error::Error(
        "Adapter '{}' resolved {} SIDs but assembled {} SID components", design.name, resolved_sids.size(), sid_index);
  }
}

}  // namespace

AdapterDesign ToAdapterDesign(const std::string& name, const CompiledAdapter& adapter) {
  // Selection-time (strict) validation: a design is only ever assembled here, when it is actually
  // selected for a run, so this is where the strict duplex arm-symmetry contract must hold.
  ValidateCompiledAdapter(name, adapter, /*check_arm_symmetry=*/true);

  // Single derivation site: build the structural resolved view (no pools loaded) and annotate every
  // id with its derived bait / UMI flank affixes. The assembler reads these back; it never re-derives.
  ResolvedAdapter resolved = ResolveStructure(name, adapter);
  Infer(resolved);

  AdapterDesign design = BuildForType(name, adapter, resolved);
  AttachStructuralSidFlanks(design, resolved);
  return design;
}

}  // namespace xoos::demux::dsl
