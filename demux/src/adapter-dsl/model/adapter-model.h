#pragma once

#include <optional>
#include <string>
#include <utility>
#include <variant>

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::demux::dsl {

/**
 * @brief The single durable model for a compiled adapter (see design docs "Two forms, one model").
 *
 * There are two JSON representations but one model:
 *   - the human-authored form (`type` + `dsl` string + `adapters` map), and
 *   - the compiled form (`type` + `read_structure` token array + `adapters` map).
 * Both compile down to the structures in this header. The compiled form is the durable artifact
 * demux consumes; the DSL is only an ergonomic authoring layer.
 *
 * Design rules encoded here (verbatim from the docs):
 *   - Case = token class: an UPPERCASE token in read_structure is a fixed literal; a lowercase
 *     token is an id that resolves in `adapters` (or the reserved `insert` keyword).
 *   - +/- = orientation: a trailing +/- on an id declares strand. + and no marker are equivalent.
 *   - Matching base names link: `sid+` and `sid-` share base `sid` and must match (RC-validated).
 *   - Nested arrays = the duplex fold: each `(...)` arm compiles to one nested token array; the
 *     fixed literal between the two arrays is the loop, recovered by position.
 */

/**
 * @brief Strand orientation of an id token relative to the FASTA it links to.
 * `+` (or no marker) is forward; `-` is reverse complement.
 */
enum class Orientation {
  kForward,
  kReverse,
};

/**
 * @brief Element type of an `adapters` entry.
 *
 * `sid` and `umi` entries describe FASTA-backed candidate pools. `fixed` entries are optional
 * side-cars for uppercase fixed literals in the read structure; they are validated at config time
 * and ignored by runtime matcher construction.
 */
enum class EntryType {
  kSid,
  kUmi,
  kFixed,
};

/**
 * @brief The adapter design family a compiled adapter belongs to.
 *
 * This is the typed form of the design `type` string carried in every config entry. The type is
 * intrinsic to a design -- it selects the code path (matcher family, metrics, formatting) the
 * runtime uses -- so it is captured as a closed enum rather than a free-form string, and unknown
 * values are rejected at parse time (see @ref ParseDesignType).
 *
 * The variants mirror the runtime @ref xoos::demux::AdapterType one-to-one; the assembler maps between
 * them. The `duplex_umi` <- `duplex_stem` <- `duplex` specialization is expressed in the runtime
 * demux-object hierarchy, not here: at the model layer these are three sibling design families that
 * each select their own compile path.
 */
enum class DesignType {
  kDuplex,
  kDuplexStem,
  kDuplexUmi,
  kYsu,
  kYs,
  kSimplex,
  kSimplex10x,
};

/// @brief A left/right wiggle window (v1 legacy search knob).
struct Wiggle {
  s64 left{0};
  s64 right{0};

  friend bool operator==(const Wiggle& a, const Wiggle& b) = default;
};

/**
 * @brief A fixed literal sequence token (UPPERCASE ACGT). This is the data itself.
 *
 * A fixed token in `read_structure` is always a bare string. It may have an optional `type: fixed`
 * side-car entry in `adapters`, but no lookup-table/search parameters live on fixed components.
 * Downstream assembly turns this token into the fixed sequence representation expected by the
 * selected trimmer.
 */
struct FixedToken {
  std::string sequence;

  friend bool operator==(const FixedToken& a, const FixedToken& b) = default;
};

/**
 * @brief A named id token (sid / umi) with an explicit orientation.
 * `base` is the id without the +/- suffix and is the key into the `adapters` map.
 */
struct IdToken {
  std::string base;
  Orientation orientation{Orientation::kForward};

  friend bool operator==(const IdToken& a, const IdToken& b) = default;
};

/**
 * @brief The reserved insert keyword. Self-describing: no id, no `adapters` entry.
 * Orientation is optional (a bare `[ins]` carries none; `[ins+]` / `[ins-]` set it).
 */
struct InsertToken {
  std::optional<Orientation> orientation;

  friend bool operator==(const InsertToken& a, const InsertToken& b) = default;
};

/// @brief A leaf token in the read structure: fixed literal, id, or the insert keyword.
using LeafToken = std::variant<FixedToken, IdToken, InsertToken>;

/**
 * @brief A duplex arm: one `(...)` group, compiled to a nested array of leaf tokens.
 * Arms never nest; an arm contains only leaf tokens.
 */
struct ArmToken {
  vec<LeafToken> tokens;

  friend bool operator==(const ArmToken& a, const ArmToken& b) = default;
};

/**
 * @brief A top-level token in the read structure: a leaf token or a duplex arm (nested array).
 * The fixed literal that appears (as a FixedToken) between two ArmTokens is the duplex loop,
 * identified by position, not by a type tag.
 */
using Token = std::variant<FixedToken, IdToken, InsertToken, ArmToken>;

/**
 * @brief Per-element configuration for a named id, keyed by base name in the `adapters` map.
 *
 * For sid/umi, `sequence` is usually a FASTA path and the remaining fields are v1 legacy matcher
 * knobs carried for parity with the shipping matcher. Sample-sheet-only SID entries can omit
 * `sequence`. For fixed, `sequence` is normalized to the fixed literal key and every matcher knob
 * must remain unset.
 */
struct AdapterEntry {
  EntryType type{EntryType::kSid};
  /// FASTA path to the candidate pool for sid/umi; empty for sample-sheet-only SID entries; literal bases for fixed.
  std::string sequence;

  std::optional<s64> max_edit_distance;
  /// `wiggle+`
  std::optional<Wiggle> wiggle_forward;
  /// `wiggle-`
  std::optional<Wiggle> wiggle_reverse;
  /// `bait_length+`
  std::optional<s64> bait_length_forward;
  /// `bait_length-`
  std::optional<s64> bait_length_reverse;
  /// umi extension
  std::optional<s64> additional_prefix_length;
  /// umi extension
  std::optional<s64> additional_suffix_length;

  friend bool operator==(const AdapterEntry& a, const AdapterEntry& b) = default;
};

/**
 * @brief A single compiled adapter: a free-form design type, an ordered read structure, and
 * the map of named id entries.
 *
 * `adapters` preserves insertion order (vector of pairs) so serialization is stable and matches
 * the order ids first appear in the read structure.
 */
struct CompiledAdapter {
  /// The design family this adapter belongs to; selects the runtime code path.
  DesignType type{DesignType::kSimplex};
  vec<Token> read_structure;
  vec<std::pair<std::string, AdapterEntry>> adapters;

  /// @brief Look up an entry by base id. Returns nullptr if absent.
  const AdapterEntry* Find(std::string_view base) const;

  /// @brief True if the read structure contains at least one duplex arm.
  bool IsDuplex() const;

  friend bool operator==(const CompiledAdapter& a, const CompiledAdapter& b) = default;
};

/// @brief A named collection of compiled adapters (top-level container / bundle).
struct AdapterBundle {
  vec<std::pair<std::string, CompiledAdapter>> adapters;

  /// @brief Look up an adapter by name. Throws xoos::Error if absent.
  const CompiledAdapter& At(const std::string& name) const;

  /// @brief True if an adapter with the given name exists.
  bool Contains(const std::string& name) const;  // NOSONAR(S5536) -- exercised by tests (not in sonar.sources)

  friend bool operator==(const AdapterBundle& a, const AdapterBundle& b) = default;
};

// --- Enum <-> string helpers (used by JSON and error messages) ---

std::string_view ToString(EntryType type);
EntryType ParseEntryType(std::string_view str);

/// @brief The canonical `type` string for a design family (e.g. "duplex", "simplex-10x").
std::string_view ToString(DesignType type);

/**
 * @brief Parse a design `type` string into a DesignType.
 * @throws xoos::Error if @p str is not one of the known design families. The message lists the
 * valid values so a config typo is actionable.
 */
DesignType ParseDesignType(std::string_view str);

/// @brief The orientation suffix character for an id token ("+", "-").
char OrientationSuffix(Orientation orientation);

/// @brief True if every character is an uppercase ACGT base (a valid fixed literal).
bool IsUpperAcgt(std::string_view s);

/// @brief Regex check for a valid lowercase id base: ^[a-z][a-z0-9_]+$.
bool IsValidIdBase(std::string_view s);

/**
 * @brief Encode a leaf token as its compiled-form string (e.g. "CAACAAT", "sid+", "insert").
 * This is the single place the case=class and +/- orientation conventions are written out.
 */
std::string LeafToCompiledString(const LeafToken& leaf);

/**
 * @brief Decode a compiled-form read_structure string token into a leaf token.
 * Applies the case=class rule: UPPERCASE ACGT -> fixed; the `insert` keyword -> insert;
 * otherwise a lowercase id with an optional +/- suffix. Throws xoos::Error on a malformed token.
 */
LeafToken LeafFromCompiledString(std::string_view token);

/**
 * @brief Validate a compiled adapter, throwing xoos::Error with an actionable message on any
 * violation. Checks: case=class, id resolution, insert rules, duplex arm count + loop-by-position,
 * and shared-name linkage consistency.
 *
 * @param check_arm_symmetry When true, also enforces strict duplex arm symmetry: the two arms must
 * be exact reverse complements (mirrored fixed literals and paired ids on opposite strands). This
 * is the strict, selection-time check applied when an adapter is actually chosen for a run. The
 * bundle load path passes false so that a config may legitimately carry designs with asymmetric
 * arms (they load, but throw here if and only if they are selected).
 */
void ValidateCompiledAdapter(const std::string& name, const CompiledAdapter& adapter, bool check_arm_symmetry);

/**
 * @brief Validate a compiled adapter with strict arm symmetry enforced (the selection-time check).
 * Equivalent to ValidateCompiledAdapter(name, adapter, true).
 */
void ValidateCompiledAdapter(const std::string& name, const CompiledAdapter& adapter);

/**
 * @brief Validate every adapter in the bundle (plus name uniqueness).
 *
 * Structural validation only: arm symmetry is deferred to selection time (see
 * ValidateCompiledAdapter) so that a bundle may contain designs whose arms are not exact reverse
 * complements. Such designs load successfully and only fail when selected for a run.
 */
void ValidateBundle(const AdapterBundle& bundle);

}  // namespace xoos::demux::dsl
