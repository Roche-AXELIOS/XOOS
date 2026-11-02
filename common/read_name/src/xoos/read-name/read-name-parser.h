#pragma once

#include <string_view>

#include <xoos/types/int.h>

#include "xoos/read-name/umi.h"

namespace xoos::read_name {

/**
 * Determine if the UMI format used in the qname is the legacy format, which uses '|' as a delimiter for UMIs, or the
 * newer format, which uses ':' as a delimiter for UMIs.  Both use '|' to separate the read info and sid from the UMIs.
 *
 * @param qname The query name string to check for UMI format. It is expected to contain at least one '|' character
 * separating the read info and sid from the UMIs.
 * @return true if the qname uses the legacy UMI format ('|' delimiter), false if it uses the newer UMI format (':'
 * delimiter).
 * @throw  error::Error if the qname does not contain at least one '|' character, which is expected to separate the read
 * info and sid from the UMIs.
 */
bool IsLegacyReadNameFormat(const std::string_view& qname);

/**
 * Determine whether the read name format carries UMI fields, independent of whether a UMI value is
 * actually present on this particular read. This is a structural (format-level) check that is stable
 * across reads in the same dataset: a UMI-prep dataset emits the UMI fields (with `*` placeholders
 * for unclassified UMIs) even on partial reads, while a non-UMI dataset omits them entirely.
 *
 * - Legacy format (`<read_info>|<sid>|<umi5p>|<umi3p>`) always carries UMI fields.
 * - Newer format carries UMI fields when the section after the last '|' contains a ':'
 *   (`<sid>:<umi5p>:<umi3p>`); a bare `<sid>` after the last '|' indicates no UMI fields.
 *
 * @param qname The query name string to inspect.
 * @return true if the read name format carries UMI fields, false otherwise (including when no '|' is present).
 */
bool HasUmiFields(const std::string_view& qname) noexcept;

struct AdapterInfo {
  bool has_5p_sid{};
  bool has_3p_sid{};
  bool has_5p_umi{};
  bool has_3p_umi{};
  Umi umi5p;
  Umi umi3p;
};

constexpr u8 kAdapterFlagHas3pUmi = 1;
constexpr u8 kAdapterFlagHas5pUmi = 2;
constexpr u8 kAdapterFlagHas3pSid = 4;
constexpr u8 kAdapterFlagHas5pSid = 8;

/**
 * @brief Parses a query name (qname) string to extract two UMI (Unique Molecular Identifier) values
 * and the bitfield representing the presence of SIDs and UMIs.
 *
 * This function supports two formats:
 * - Legacy: <read_info>|<sid>|<umi5p>|<umi3p> (only supports UMI parsing)
 * - Newer:  <read_info>:<bitfield>|<sid>:<umi5p>:<umi3p> (supports both UMI parsing and presence bitfield parsing)
 *
 * This function handles the following scenarios:
 * 1. UMIs Present: Extracts the 5' and 3' UMI sequences.
 * 2. UMIs Missing: If either UMI is missing (represented by `kNoUmi`), it is returned as `std::nullopt`.
 * 3. UMIs Omitted: If the dataset does not support UMIs and the fields are omitted entirely
 * (e.g., `<read_info>|<sid>`), this function will throw an error because the expected delimiters
 * are not found and the function expects to only be called for read names with UMIs.
 *
 * @param qname The input string containing the query name with embedded UMIs, separated by '|'.
 * @param parse_umi A boolean indicating whether to parse UMIs from the read name. If true, the function will attempt to
 *                    parse UMIs from the read name and return them in the output. If false, the function will skip UMI
 *                    parsing and return `std::nullopt` for both UMIs in the output.
 * @return A struct containing the parsed UMI values and the presence information for SIDs and UMIs.
 * @throws error::Error If read name is ill-formed or does not follow one of the two expected formats.
 */
AdapterInfo ParseReadName(const std::string_view& qname, bool parse_umi);

/**
 * @brief A no-throw wrapper around `ParseReadName` that returns default values in case of parsing errors.
 *
 * This function calls `ParseReadName` to extract UMI information from the query name. If `ParseReadName` throws an
 * exception due to parsing errors (e.g., ill-formed read name, missing delimiters), this function catches the exception
 * and returns a default `AdapterInfo` object with the following values:
 * - `has_5p_sid` and `has_3p_sid` set to true, indicating the presence of SIDs.
 * - `has_5p_umi` and `has_3p_umi` set to false, indicating the absence of UMIs.
 * - `umi5p` and `umi3p` set to `std::nullopt`, indicating no UMI values.
 *
 * @param qname The input string containing the query name with embedded UMIs, separated by '|'.
 * @param parse_umi A boolean indicating whether to parse UMIs from the read name. If true, the function will attempt to
 *                    parse UMIs from the read name and return them in the output. If false, the function will skip UMI
 *                    parsing and return `std::nullopt` for both UMIs in the output.
 * @return An `AdapterInfo` object containing UMI presence information and UMI values (if parsed successfully).
 */
AdapterInfo ParseReadNameNoExcept(const std::string_view& qname, bool parse_umi) noexcept;

}  // namespace xoos::read_name
