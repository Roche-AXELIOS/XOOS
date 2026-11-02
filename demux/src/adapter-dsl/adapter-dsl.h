#pragma once

/**
 * @file adapter-dsl.h
 * @brief Public facade for the DSL adapter configuration: DSL authoring, a compiled JSON
 * artifact, and runtime resolution.
 *
 * Two forms, one model (see the demux custom-adapter design docs):
 *
 *   authored JSON  --FromAuthoredJson-->  model  --CompiledToJson-->  compiled JSON
 *   { type, dsl,                          (AdapterBundle)             { type,
 *     adapters }                                                        read_structure,
 *                                                                       adapters }
 *
 *   compiled JSON  --CompiledFromJson-->  model  --Resolve-->  runtime elements
 *                                                              (ResolvedAdapter, FASTA pools)
 *
 * The DSL is the ergonomic authoring layer; the compiled JSON is the durable artifact demux
 * consumes; resolution loads and orients the candidate pools for the matcher. All entry points
 * validate loudly (see ValidateBundle).
 */

// This is a facade header: it deliberately re-exports the per-layer public headers so callers can
// include a single file. The IWYU export pragmas tell include-cleaner these re-exports are
// intentional, not unused includes.
// IWYU pragma: begin_exports
#include "adapter-dsl/authored/authored-json.h"        // FromAuthoredJson, FromAuthoredJsonFile
#include "adapter-dsl/dsl/dsl-compiler.h"              // CompileDsl, CompileDslToBundle
#include "adapter-dsl/model/adapter-model.h"           // model + validation
#include "adapter-dsl/model/compiled-json.h"           // CompiledToJson, CompiledFromJson(File)
#include "adapter-dsl/parser/read-structure-parser.h"  // ReadStructureParser, ParserRegistry
#include "adapter-dsl/resolve/resolved-adapter.h"      // Resolve, ResolvedAdapter
// IWYU pragma: end_exports
