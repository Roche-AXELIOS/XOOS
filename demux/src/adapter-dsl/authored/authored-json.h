#pragma once

#include <filesystem>  // NOLINT(build/c++17)
#include <string>

#include <xoos/vfs/vfs.h>

#include "adapter-dsl/model/adapter-model.h"

namespace fs = std::filesystem;

namespace xoos::demux::dsl {

/**
 * @brief Load the human-authored JSON form and compile it into the model.
 *
 * The authored form pairs a single DSL string with the linkage/params it cannot express inline
 * (FASTA paths and v1 legacy tuning knobs):
 *
 *   {
 *     "YSU": {
 *       "type": "YSU",
 *       "dsl": "CAACAAGAGTCTTTT[sid:sid-]GCG[umi:umi1+]GAGTT[ins]AACTC[umi:umi2-]CGC[sid:sid+]...",
 *       "adapters": {
 *         "sid":  { "type": "sid", "sequence": "sids/hpsid.fa", "max_edit_distance": 2 },
 *         "umi1": { "type": "umi", "sequence": "ysu/umi-5p.fa" },
 *         "umi2": { "type": "umi", "sequence": "ysu/umi-5p.fa" }
 *       }
 *     }
 *   }
 *
 * Compilation: the DSL yields the read_structure and a bare adapters skeleton (ids + types); the
 * authored `adapters` block then overlays sequence and legacy params onto that skeleton. The
 * result is a fully-populated, validated compiled adapter identical in shape to what
 * compiled-json produces, so `CompiledToJson(FromAuthoredJson(...))` is the DSL->compiled step.
 */

/// @brief Compile an authored JSON bundle text into the model. Throws xoos::Error on any error.
AdapterBundle FromAuthoredJson(const std::string& json_text);

/// @brief Compile an authored JSON bundle file into the model. NOSONAR(S5536) exercised by tests only.
AdapterBundle FromAuthoredJsonFile(const fs::path& path);  // NOSONAR(S5536)

/**
 * @brief Compile an authored JSON bundle read from inside an adapter design bundle (ZIP or directory).
 *
 * @param vfs  An open virtual filesystem rooted at the adapter design bundle.
 * @param path Path of the config file within the bundle (e.g. "adapter-config.json5").
 * FASTA pools referenced by the config remain relative to the bundle and are resolved downstream.
 */
AdapterBundle FromAuthoredJsonBundle(const vfs::VirtualFilesystemPtr& vfs, const fs::path& path);

}  // namespace xoos::demux::dsl
