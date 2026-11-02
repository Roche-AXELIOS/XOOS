#pragma once

#include <optional>
#include <string>
#include <string_view>

#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

namespace xoos::io {

/**
 * @name Canonical @RG string values for SBX BAMs
 *
 * The wire values written into @RG tags, exposed as a single source of truth so
 * production code, downstream tools, and tests never re-spell them.
 *
 * The run_type values (kRunType*) are not modelled as an enum on purpose. The
 * I2X_* families use versioned names, so LibraryTypeFromRunType matches their
 * stable prefixes case-insensitively.
 * @{
 */

/// @RG PL (platform) value stamped by the xoos pipeline for all SBX data.
inline constexpr std::string_view kSbxPlatform = "SBX";
/// @RG PM (platform model) value: the SBX instrument model on newer BAMs. Older
/// instrument BAMs (ppa2bam) instead carry the run_type in PM; see DetectLibraryType.
inline constexpr std::string_view kSbxPlatformModel = "AXELIOS_1";

// Duplex run_types. SBX-FAST is its own representative; the rest map to kSbxD.
inline constexpr std::string_view kRunTypeSbxFast = "SBX-FAST";
inline constexpr std::string_view kRunTypeSbxD = "SBX-D";
inline constexpr std::string_view kRunTypeSbxDm = "SBX-DM";
inline constexpr std::string_view kRunTypeI2xI7HpsidDuplexV0 = "I2X_i7_HPSID_duplex_v0";
inline constexpr std::string_view kRunTypeI2xHpsidPrefix = "I2X_HPSID";
inline constexpr std::string_view kRunTypeI2xHpuPrefix = "I2X_HPU";

// Simplex run_types, all mapping to kSbxS.
inline constexpr std::string_view kRunTypeYs = "YS";
inline constexpr std::string_view kRunTypeYsu = "YSU";
inline constexpr std::string_view kRunTypeYsNew = "YS-NEW";

/// @}

/**
 * @brief Library type of an SBX BAM file, as three representative classes.
 *
 * All SBX data shares the same sequencing chemistry, so what varies -- and what
 * this enum captures -- is the library prep protocol / run type. The run_type is
 * stamped into an @RG tag: the xoos pipeline sets LB to run_type verbatim (see
 * the samplesheet run_type enum and subworkflows/local/align/nextflow.config),
 * while raw instrument BAMs from ppa2bam set PM instead (e.g. "PM:SBX-D", with
 * LB holding the flow cell/run id). Downstream tools only need to know which of
 * three behavioural classes the data falls into. The many run_type values are
 * therefore collapsed onto these representatives (see LibraryTypeFromRunType):
 *
 *   - kSbxFast <- SBX-FAST (its own class: duplex, but distinct handling).
 *   - kSbxD    <- every other duplex run_type: SBX-D, SBX-DM,
 *                 I2X_i7_HPSID_duplex_v0, and designs whose names begin with
 *                 I2X_HPSID or I2X_HPU.
 *   - kSbxS    <- every simplex run_type: YS, YSU, YS-NEW.
 */
enum class LibraryType : u8 {
  /// Simplex library (single read from one strand). Represents YS, YSU, YS-NEW.
  kSbxS,
  /// Duplex library. Represents SBX-D, SBX-DM, and the I2X_* designs.
  kSbxD,
  /// SBX-FAST duplex library, kept distinct from the other duplex designs.
  kSbxFast,
};

/// Number of reads sampled from the start of the BAM for the content-based
/// detection fallbacks (YC-tag scan and base-quality heuristic).
inline constexpr u32 kLibraryTypeReadSampleSize = 10'000;

/// Fraction of sampled bases with BQ > 30 above which a dataset is considered
/// not simplex. Duplex concordant bases carry high legacy base quality.
inline constexpr f64 kLibraryTypeHighBqFraction = 0.05;

/// Base-quality threshold used by the high-BQ heuristic.
inline constexpr u8 kLibraryTypeHighBqThreshold = 30;

/**
 * @brief Map a pipeline run_type string to its representative LibraryType.
 *
 * Matching is case-insensitive. The accepted run_type values are the samplesheet
 * schema enum (subworkflows align/nextflow.config stamps run_type into @RG LB
 * verbatim), collapsed onto the three representatives:
 *   - "SBX-FAST"                                     -> kSbxFast
 *   - "SBX-D", "SBX-DM", "I2X_i7_HPSID_duplex_v0", prefixes "I2X_HPSID",
 *     "I2X_HPU"                                        -> kSbxD (other duplex)
 *   - "YS", "YSU", "YS-NEW"                          -> kSbxS (simplex)
 *
 * Values of "unknown" or empty, and any unrecognised value, return std::nullopt.
 * Prefix matching makes versioned I2X_* designs forward-compatible and allows
 * their prefix component to be recognized after splitting a structured field.
 *
 * @param run_type The run_type / @RG LB (or legacy PM) tag value.
 * @return The representative LibraryType, or std::nullopt when unrecognised.
 */
std::optional<LibraryType> LibraryTypeFromRunType(const std::string& run_type);

/**
 * @brief Detect the library type of a BAM file.
 *
 * Resolution order (stops at the first match; the returned type is logged in
 * every case):
 *   1. @RG LB component 1 -- adapter class from the structured library value.
 *                     Logged at INFO.
 *   2. Full @RG LB tag -- run_type as stamped by older xoos pipelines. Logged at
 *                     INFO.
 *   3. @RG PU component 4 -- adapter class from the structured platform unit.
 *                     Logged at INFO.
 *   4. Full @RG PM tag -- run_type as stamped by the instrument. ppa2bam (the
 *                     on-instrument tool) writes the run_type into PM (e.g.
 *                     "PM:SBX-D") and uses LB for the flow cell/run id, so raw
 *                     instrument BAMs only resolve here; some legacy
 *                     germline_wgs BAMs likewise carry the run_type in PM.
 *                     Logged at WARNING (PU and LB did not resolve the run_type).
 *   5. YC-tag scan -- the first ::kLibraryTypeReadSampleSize reads are examined;
 *                     if any carry a YC tag encoding a duplex read the file is
 *                     classified as duplex. Logged at WARNING (inferred from
 *                     read content, not the header).
 *   6. Simplex BQ heuristic -- if none of the sampled reads carry a YC tag and
 *                     either one distinct base quality is observed, or two to
 *                     four distinct qualities are observed and each is one of
 *                     Q10, Q14, Q15, and Q20, the file is classified as simplex.
 *                     Logged at WARNING.
 *   7. Duplex BQ heuristic -- if more than ::kLibraryTypeHighBqFraction of sampled bases
 *                     exceed ::kLibraryTypeHighBqThreshold the file is not simplex,
 *                     so it is classified as duplex. Logged at WARNING.
 *
 * The header modes (1 through 4) map the extracted or full value onto a
 * representative via ::LibraryTypeFromRunType. Structured LB has priority over
 * full LB, which has priority over structured PU and full PM. The duplex content
 * fallbacks (5 and 7) cannot
 * distinguish SBX-D from SBX-FAST, so they return ::LibraryType::kSbxD (the
 * representative for non-FAST duplex). The simplex content fallback (6) returns
 * ::LibraryType::kSbxS only when no sampled read carries any YC tag.
 *
 * Tag values of "unknown" (case-insensitive) or empty are treated as non-matches
 * so detection falls through to the next mode. Missing tags never throw.
 *
 * @param bam_path Path to the BAM/CRAM file.
 * @return The detected LibraryType, or std::nullopt when it cannot be resolved
 *         (logged at WARNING).
 */
std::optional<LibraryType> DetectLibraryType(const fs::path& bam_path);

/**
 * @brief Whether a BAM file's library type is duplex.
 *
 * @param bam_path Path to the BAM/CRAM file.
 * @return true for duplex library types (kSbxD, kSbxFast), false for simplex
 *         (kSbxS), std::nullopt when the library type cannot be detected.
 */
std::optional<bool> IsDuplex(const fs::path& bam_path);

/// True for duplex library types (kSbxD, kSbxFast), false for simplex (kSbxS).
bool IsDuplex(LibraryType library_type);

/**
 * @brief Cross-check a header-declared library type against the read content.
 *
 * Detects contradictions between what the header claims and what the reads show,
 * over the first ::kLibraryTypeReadSampleSize reads:
 *   - declared duplex (kSbxD / kSbxFast) but the reads carry no YC tags;
 *   - declared simplex (kSbxS) but the reads carry duplex-encoded YC tags, or
 *     more than ::kLibraryTypeHighBqFraction of bases exceed
 *     ::kLibraryTypeHighBqThreshold.
 *
 * This reports the mismatch; it does not override the declared library type.
 * An error is logged for each contradiction.
 *
 * @param bam_path Path to the BAM/CRAM file.
 * @param declared The library type declared by the header (e.g. from ::DetectLibraryType).
 * @return std::nullopt when the declared library type is consistent with the read
 *         content; otherwise a message describing the inconsistency.
 */
std::optional<std::string> CheckLibraryTypeConsistency(const fs::path& bam_path, LibraryType declared);

}  // namespace xoos::io
