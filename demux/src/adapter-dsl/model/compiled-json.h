#pragma once

#include <filesystem>  // NOLINT(build/c++17)
#include <string>

#include <xoos/types/int.h>

#include "adapter-dsl/model/adapter-model.h"

namespace fs = std::filesystem;

namespace xoos::demux::dsl {

/**
 * @brief (De)serialization for the compiled JSON form.
 *
 * The compiled form is the durable artifact demux consumes. Its shape (from the design docs):
 *
 *   {
 *     "SBX-D": {
 *       "type": "duplex",
 *       "read_structure": [
 *         "CAACAAT",
 *         ["GCACCTGACGTACGCT", "insert+", "A", "sid+"],
 *         "ATGTATG",
 *         ["sid-", "T", "insert-", "AGCGTACGTCAGGTGC"],
 *         "TGAAGAGC"
 *       ],
 *       "adapters": { "sid": { "type": "sid", "sequence": "sids/hpsid.fa" } }
 *     }
 *   }
 *
 * A read_structure token is a plain string when it is a leaf (fixed literal / `sid+` id / the
 * `insert` keyword) and a nested array when it is a duplex arm. Case tells fixed from id; the
 * trailing +/- gives orientation. `adapters` re-materializes the `type` on each entry.
 */

/**
 * @brief Serialize a bundle to canonical compiled JSON text (pretty-printed).
 * @brief Serialize @p bundle to compiled-form JSON with the given @p indent width.
 */
std::string CompiledToJson(const AdapterBundle& bundle, s32 indent);

/// @brief Serialize @p bundle to compiled-form JSON with the default indent width (2 spaces).
std::string CompiledToJson(const AdapterBundle& bundle);  // NOSONAR(S5536) -- exercised by tests (not in sonar.sources)

/// @brief Parse a bundle from compiled JSON text. Throws xoos::Error on malformed input.
AdapterBundle CompiledFromJson(const std::string& json_text);

/// @brief Parse a bundle from a compiled JSON file path.
AdapterBundle CompiledFromJsonFile(  // NOSONAR(S5536) -- exercised by tests (not in sonar.sources)
    const fs::path& path);

}  // namespace xoos::demux::dsl
