#pragma once

#include <string>
#include <variant>

#include <nlohmann/json_fwd.hpp>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "sequence/sequence-enums.h"

namespace xoos::demux {
/**
 * @brief The class of a variable adapter component that builds a sequence matcher.
 *
 * Fixed literals do not have a matcher kind: they are represented by @ref FixedComponentDefinition.
 */
enum class MatcherKind {
  kUmi,
  kSid,
};

/**
 * @brief Positional search window for a variable matcher around its expected location.
 */
struct SearchWindow {
  s32 max_wiggle_left{0};
  s32 max_wiggle_right{0};

  friend bool operator==(const SearchWindow& a, const SearchWindow& b) = default;
};

/**
 * @brief Constant adapter literal in molecule read order.
 *
 * Fixed components are consumed directly by trimmers; no LUT or matcher settings are stored here.
 */
struct FixedComponentDefinition {
  std::string sequence;

  friend bool operator==(const FixedComponentDefinition& a, const FixedComponentDefinition& b) = default;
};

/**
 * @brief Variable adapter component that is materialized as a @ref SeqMatcher.
 *
 * `prefix` / `suffix` are literal affixes. When a UMI or SID needs a flank extension, the DSL
 * inference layer derives it positionally and writes the literal value here.
 */
struct MatcherComponentDefinition {
  MatcherKind kind{};
  /// @brief FASTA file containing the candidate barcode pool. Empty for sample-sheet-only SID components.
  fs::path sequences;
  /// @brief Maximum edit distance expanded into the matcher LUT.
  u32 max_edit_distance{0};
  /// @brief Transform applied after prefix/suffix affixes and before matcher construction.
  SeqTransform transform{SeqTransform::kNone};
  /// @brief Bases prepended to each FASTA sequence before transform and matcher construction.
  std::string prefix;
  /// @brief Bases appended to each FASTA sequence before transform and matcher construction.
  std::string suffix;

  /**
   * @brief Runtime-only DSL-derived sequence of this SID's preceding element when that neighbor is a
   * fixed literal (else empty), in read orientation. Regenerated during assembly and not serialized.
   *
   * TODO: Reference the read structure by index rather than storing these raw flank copies. Once
   * AdapterDesign carries its ordered elements (including the insert) -- mirroring ResolvedId's
   * preceding_element / following_element -- a matcher can point at its neighbor and derive the flank
   * on demand, so structural_prefix / structural_suffix become computed rather stored here
   */
  std::string structural_prefix;
  /// @brief As @ref structural_prefix, for the following element (immediately after this SID).
  std::string structural_suffix;

  /// @brief Positional search window around the expected component location.
  SearchWindow search{};
  /// @brief Optional bait seed for a barcode (sid): the adjacent fixed trimmed to the configured
  /// bait length, in read order. Empty when the component has no bait. Attached to the component
  /// itself so it is addressed by position, not by a 5'/3' side tag.
  std::string bait;

  friend bool operator==(const MatcherComponentDefinition& a,
                         const MatcherComponentDefinition& b) = default;  // NOSONAR(S1067)
};

using AdapterComponentDefinition = std::variant<FixedComponentDefinition, MatcherComponentDefinition>;

std::string FormatJson(const AdapterComponentDefinition& component);

enum class AdapterType {
  kDuplex,
  kDuplexUMI,
  kDuplexStem,
  kYsu,
  kYs,
  kSimplex,
  kSimplex10x
};

struct AdapterDesign {
  std::string name;
  AdapterType type{};
  /// @brief All adapter components in exact DSL read_structure order (5' end to 3' end of the
  /// molecule, duplex arms concatenated with the loop as the natural boundary). There is no 5'/3'
  /// split: the trimming algorithm addresses every component by its position in this vector and
  /// decides how to use each. Orientation is baked into each component's transform by the assembler.
  vec<AdapterComponentDefinition> components;
};

/**
 * Load the named adapter design from the given bundle.
 *
 * The bundle (ZIP or directory) must contain a DSL adapter config ("adapter-config.json5"). The named
 * design is compiled and assembled into an AdapterDesign; FASTA pools referenced by the config are resolved
 * relative to the bundle downstream. Empty SID pool paths must be supplied by a sample sheet override.
 */
AdapterDesign LoadAdapterDesign(const fs::path& adapter_design_bundle, const std::string& adapter_design_name);

void to_json(nlohmann::json& j, const MatcherKind& e);                   // NOSONAR(S100)
void from_json(const nlohmann::json& j, MatcherKind& e);                 // NOSONAR(S100)
void to_json(nlohmann::json& j, const SeqTransform& e);                  // NOSONAR(S100)
void from_json(const nlohmann::json& j, SeqTransform& e);                // NOSONAR(S100)
void to_json(nlohmann::json& j, const AdapterType& e);                   // NOSONAR(S100)
void from_json(const nlohmann::json& j, AdapterType& e);                 // NOSONAR(S100)
void to_json(nlohmann::json& j, const SearchWindow& d);                  // NOSONAR(S100)
void from_json(const nlohmann::json& j, SearchWindow& d);                // NOSONAR(S100)
void to_json(nlohmann::json& j, const FixedComponentDefinition& d);      // NOSONAR(S100)
void from_json(const nlohmann::json& j, FixedComponentDefinition& d);    // NOSONAR(S100)
void to_json(nlohmann::json& j, const MatcherComponentDefinition& d);    // NOSONAR(S100)
void from_json(const nlohmann::json& j, MatcherComponentDefinition& d);  // NOSONAR(S100)
void to_json(nlohmann::json& j, const AdapterComponentDefinition& d);    // NOSONAR(S100)
void from_json(const nlohmann::json& j, AdapterComponentDefinition& d);  // NOSONAR(S100)
void to_json(nlohmann::json& j, const AdapterDesign& d);                 // NOSONAR(S100)
void from_json(const nlohmann::json& j, AdapterDesign& d);               // NOSONAR(S100)

}  // namespace xoos::demux
