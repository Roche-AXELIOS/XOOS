#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "adapter-dsl/model/adapter-model.h"

namespace xoos::demux::dsl {

/**
 * @file read-structure-parser.h
 * @brief The pluggable seam for authoring syntaxes.
 *
 * The model, authored-JSON, compiled-JSON, and resolve layers are all parser-agnostic: they only
 * ever see a @ref CompiledAdapter. The *only* thing tied to a concrete syntax is the step that
 * turns an authored string into that model. This header makes that step swappable.
 *
 * A @ref ReadStructureParser turns one authored-syntax string (today: the DSL) into a
 * @ref CompiledAdapter carrying the read_structure and a bare adapters skeleton (ids + types). The
 * authored-JSON layer then overlays FASTA paths and params. Parsers are registered by a syntax
 * name in a @ref ParserRegistry, and the authored form selects one via an optional `syntax` field
 * (defaulting to @ref kDefaultSyntax). Adding a new frontend (PEG, lexy, S-expr, ...) is a matter
 * of implementing this interface and registering it in one place (builtin-parsers.cpp); no call
 * site changes.
 */

/// @brief The default authoring syntax when an authored adapter does not name one.
constexpr std::string_view kDefaultSyntax = "dsl";

/// @brief JSON key selecting the authoring syntax in the authored form (optional).
constexpr std::string_view kSyntaxKey = "syntax";

/**
 * @brief Turns one authored-syntax string into the compiled model (read_structure + skeleton).
 *
 * Implementations must be stateless and thread-safe: a single instance is shared via the registry.
 * They produce only the read_structure and a bare adapters map (one empty-sequence entry per
 * distinct id, with its type set); the authored-JSON layer overlays the rest.
 */
class ReadStructureParser {
 public:
  ReadStructureParser() = default;
  ReadStructureParser(const ReadStructureParser&) = delete;
  ReadStructureParser(ReadStructureParser&&) = delete;
  ReadStructureParser& operator=(const ReadStructureParser&) = delete;
  ReadStructureParser& operator=(ReadStructureParser&&) = delete;
  virtual ~ReadStructureParser() = default;

  /// @brief The syntax name this parser handles (e.g. "dsl"). Used as the registry key.
  virtual std::string_view Syntax() const = 0;

  /// @brief Parse @p text into a compiled adapter. Throws xoos::Error on any syntax error.
  virtual CompiledAdapter Parse(std::string_view text) const = 0;
};

/**
 * @brief Registry mapping a syntax name to its parser.
 *
 * The process-wide @ref Instance() is a shared, read-only registry with the built-in parsers
 * already registered, so callers can `Instance().Get("dsl")` without any explicit setup. Because
 * that shared instance is const, @ref Register cannot mutate it; callers that need a custom set of
 * parsers (e.g. tests, or overriding a built-in) construct their own @ref ParserRegistry and
 * register onto it. Registration is idempotent-by-replacement: registering an existing name
 * overrides it.
 */
class ParserRegistry {
 public:
  ParserRegistry() = default;

  /// @brief The shared read-only registry, with built-in parsers already registered.
  static const ParserRegistry& Instance();

  /// @brief Register (or replace) the parser for its @ref ReadStructureParser::Syntax name.
  void Register(std::shared_ptr<ReadStructureParser> parser);

  /// @brief Look up the parser for @p syntax. Throws xoos::Error (listing known syntaxes) if absent.
  const ReadStructureParser& Get(std::string_view syntax) const;

 private:
  std::vector<std::pair<std::string, std::shared_ptr<ReadStructureParser>>> _parsers;
};

/// @brief Convenience: parse @p text with the parser registered for @p syntax.
CompiledAdapter ParseReadStructure(std::string_view syntax, std::string_view text);

/**
 * @brief Register every parser shipped with demux into @p registry.
 *
 * Defined in builtin-parsers.cpp and called lazily by ParserRegistry::Instance(). Kept as a free
 * function so the registry has no compile-time dependency on any concrete parser (the seam).
 */
void RegisterBuiltinParsers(ParserRegistry& registry);

}  // namespace xoos::demux::dsl
