#pragma once

#include <filesystem>
#include <string>

#include <xoos/types/vec.h>

#include "core/dual-sid-mode.h"
#include "sequence/matcher/match-info.h"

namespace xoos::demux {

struct AdapterDesign;

namespace fs = std::filesystem;

constexpr auto kSampleNameColumnHeader = "sample_name";
constexpr auto kSampleSidColumnHeader = "sample_sid";

/// @brief One normalized sample-sheet row in 5p/3p SID order.
struct SampleSheetEntry {
  std::string sid_5p;
  std::string sid_3p;
  std::string name;
};

/// @brief Normalized sample-sheet pools used by adapter matchers and metrics.
struct SampleSheetSidPools {
  /// @brief Sample-facing pool. Sequence is normalized as "sid_5p-sid_3p".
  SidPool sid_pool;
  /// @brief Unique 5p matcher pool.
  BarcodePool sid_5p_pool;
  /// @brief Unique 3p matcher pool.
  BarcodePool sid_3p_pool;
  /// @brief Maps unique per-side matcher IDs back to sample row IDs.
  DualSidMatchingPolicy matching_policy;
  /// @brief True when at least one row carries different 5p and 3p SID sequences.
  bool has_dual_sids{false};
};

using SampleSheetEntries = vec<SampleSheetEntry>;

/**
 * @brief Read a sample sheet (CSV or TSV) and return normalized entries in row order.
 *
 * The preferred format is columns "sample_sid" and "sample_name". "sample_sid" can be either a
 * single SID or a dual SID written as "sid_5p-sid_3p"; single SIDs normalize to "sid-sid".
 *
 * @param path Path to the sample sheet file.
 * @return Normalized entries preserving CSV/TSV row order.
 * @throws error::Error On duplicate SID pairs, duplicate sample names, or reserved name usage.
 */
SampleSheetEntries ReadSampleSheet(const fs::path& path);

/**
 * @brief Loads a sample sheet and converts it into SID pools with contiguous 0-based IDs.
 *
 * Each row becomes one sample-facing SID pair. The matcher pools contain unique 5p and 3p side
 * sequences, and @ref DualSidMatchingPolicy maps those side IDs back to sample row IDs. Partial
 * assignment eligibility is derived from per-end SID collisions and @p dual_sid_mode.
 *
 * @param sample_sheet Path to the sample sheet CSV/TSV file.
 * @param dual_sid_mode Runtime override for dual-SID single-ended assignment.
 * @return SID pools and matching policy. @ref SampleSheetSidPools::sid_pool preserves row order.
 */
SampleSheetSidPools LoadSampleSheetSidPools(const fs::path& sample_sheet, DualSidMode dual_sid_mode);

/// @brief Backward-compatible helper returning a matcher-ready single-SID pool.
BarcodePool LoadSampleSheet(const fs::path& sample_sheet);

/**
 * @brief Normalize a sample-sheet SID pool using canonical sequences and structure-derived legacy forms.
 *
 * Canonical SIDs are preserved as supplied. The complete pool is converted to canonical orientation only when every
 * SID is the reverse complement of a known canonical SID, or when every SID is a uniquely recognized legacy extension
 * with the same orientation and structural-flank placement. Legacy flanks must match the selected adapter design.
 * Mixed or unrecognized pools are preserved as supplied. Conversions and preserved non-canonical pools emit warnings.
 *
 * @param sample_sids Sample-sheet SID pool, modified in place only for a uniform recognized representation.
 * @param canonical_sids Known SID sequences in canonical orientation.
 * @param design Selected adapter design whose SID matcher affixes define valid structural flanks.
 */
void NormalizeSidPoolToCanonical(SidPool& sample_sids, const BarcodePool& canonical_sids, const AdapterDesign& design);

}  // namespace xoos::demux
