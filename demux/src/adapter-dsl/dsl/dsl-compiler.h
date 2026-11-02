#pragma once

#include <string>
#include <string_view>

#include "adapter-dsl/model/adapter-model.h"

namespace xoos::demux::dsl {

/**
 * @brief Compile a DSL adapter string into the compiled model (read_structure only).
 *
 * The DSL is the ergonomic authoring layer. A molecule is written once, 5' to 3'
 * (left-to-right), the way a person already thinks about it, e.g. (YSU):
 *
 *   CAACAAGAGTCTTTT[sid:sid-]GCG[umi:umi1+]GAGTT[ins]AACTC[umi:umi2-]CGC[sid:sid+]TTTGTC...
 *
 * Grammar (design docs sections 4-6):
 *   - Fixed sequence: an inline UPPERCASE [ATCG]+ run (the data itself).
 *   - Special element in square brackets: `[<type>:<id><strand>]` or the reserved `[ins]`.
 *       - type is the first 3 chars: `sid` | `umi` | `ins`.
 *       - `sid` / `umi` are followed by a colon and an id matching ^[a-z][a-z0-9_]+$.
 *       - a trailing +/- declares orientation (+ and no marker are equivalent).
 *       - `[ins]` is a reserved keyword with no id (a bare `[ins]`, or `[ins+]` / `[ins-]`).
 *   - Duplex arms: each `(...)` group compiles to one nested arm array; the fixed literal between
 *     the two arms is the loop, recovered by position.
 *
 * This produces only the read_structure and a bare adapters map (one empty-sequence entry per
 * distinct id, with type re-materialized). The authored-JSON layer overlays FASTA paths and the
 * v1 legacy params onto the adapters entries. Call ValidateCompiledAdapter afterwards for full
 * structural checks (arm symmetry, RC complementarity, etc.).
 *
 * @throws xoos::Error with an actionable message on any syntax or structural violation.
 */
CompiledAdapter CompileDsl(std::string_view dsl_text);

/// @brief Compile a DSL string into a single-adapter bundle under @p name (fully validated).
AdapterBundle CompileDslToBundle(  // NOSONAR(S5536) -- exercised by tests (not in sonar.sources)
    const std::string& name,
    std::string_view dsl_text);

}  // namespace xoos::demux::dsl
