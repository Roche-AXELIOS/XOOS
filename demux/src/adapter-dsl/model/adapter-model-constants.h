#pragma once

#include <string_view>

namespace xoos::demux::dsl {

/**
 * @brief String constants for the compiled and authored JSON forms, and for the DSL grammar.
 *
 * Two representations, one model (see the design docs "Two forms, one model"):
 *   - Authored form: `{ "<name>": { "type": ..., "dsl": "<string>", "adapters": {...} } }`
 *   - Compiled form: `{ "<name>": { "type": ..., "read_structure": [...], "adapters": {...} } }`
 * The compiled form is the durable artifact demux consumes; the DSL is the ergonomic
 * authoring layer that compiles into it.
 */

// --- Adapter-level JSON keys ---

/**
 * @brief Free-form design type carried verbatim (e.g. "duplex", "simplex", "YSU", "YS").
 * Per the tiering discussion this is inferable (nested arms => duplex) but kept in v1.
 */
constexpr std::string_view kTypeKey = "type";
/// @brief Authored form: the single DSL string that compiles into the read structure.
constexpr std::string_view kDslKey = "dsl";
/**
 * @brief Authored form (optional): names the authoring syntax; absent means the default DSL.
 * The concrete key string lives with the parser seam (kSyntaxKey in read-structure-parser.h).
 * @brief Compiled form: the ordered token array. Order is the molecule structure.
 */
constexpr std::string_view kReadStructureKey = "read_structure";
/// @brief Both forms: map keyed by base id (e.g. "sid", "umi1") to per-element config.
constexpr std::string_view kAdaptersKey = "adapters";

// --- bundle-level keys (wrapper shape) ---

/// @brief Bundle wrapper: map keyed by design name to each design object.
constexpr std::string_view kAdapterDesignsKey = "adapter_designs";
// --- adapters-entry JSON keys ---

/// @brief Entry type re-materialized in the compiled form ("sid" / "umi" / "fixed").
constexpr std::string_view kEntryTypeKey = "type";
/// @brief FASTA path for sid/umi entries; optional literal confirmation for fixed side-cars.
constexpr std::string_view kSequenceKey = "sequence";

// --- v1 legacy per-element params (see design docs, YSU worked example) ---

constexpr std::string_view kMaxEditDistanceKey = "max_edit_distance";
constexpr std::string_view kWiggleFwKey = "wiggle+";
constexpr std::string_view kWiggleRcKey = "wiggle-";
constexpr std::string_view kBaitLengthFwKey = "bait_length+";
constexpr std::string_view kBaitLengthRcKey = "bait_length-";
constexpr std::string_view kAdditionalPrefixLengthKey = "additional_prefix_length";
constexpr std::string_view kAdditionalSuffixLengthKey = "additional_suffix_length";
constexpr std::string_view kWiggleLeftKey = "left";
constexpr std::string_view kWiggleRightKey = "right";

// --- Entry type values ---

constexpr std::string_view kEntryTypeSid = "sid";
constexpr std::string_view kEntryTypeUmi = "umi";
constexpr std::string_view kEntryTypeFixed = "fixed";

// --- Design type values ---
//
// The canonical `type` strings a design may declare. These are the single source of truth for the
// design-type vocabulary: the DesignType enum parses/serializes against them, and they must stay in
// lockstep with the runtime AdapterType strings (adapter-design.cpp) so a design round-trips through
// the assembler unchanged.

constexpr std::string_view kDesignTypeDuplex = "duplex";
constexpr std::string_view kDesignTypeDuplexStem = "duplex_stem";
constexpr std::string_view kDesignTypeDuplexUmi = "duplex_umi";
constexpr std::string_view kDesignTypeYsu = "YSU";
constexpr std::string_view kDesignTypeYs = "YS";
constexpr std::string_view kDesignTypeSimplex = "simplex";
constexpr std::string_view kDesignTypeSimplex10x = "simplex-10x";

// --- Reserved compiled-form keyword ---

/// @brief Insert keyword as it appears in the compiled read_structure (self-describing, no entry).
constexpr std::string_view kInsertKeyword = "insert";

// --- DSL grammar tokens ---

/// @brief Insert keyword as written in the DSL (the 3-character form).
constexpr std::string_view kDslInsertKeyword = "ins";
/// @brief DSL bracket type prefix for a sample-ID element.
constexpr std::string_view kDslSid = "sid";
/// @brief DSL bracket type prefix for a UMI element.
constexpr std::string_view kDslUmi = "umi";

constexpr char kOrientationForward = '+';
constexpr char kOrientationReverse = '-';
constexpr char kBracketOpen = '[';
constexpr char kBracketClose = ']';
constexpr char kParenOpen = '(';
constexpr char kParenClose = ')';
constexpr char kColon = ':';

}  // namespace xoos::demux::dsl
