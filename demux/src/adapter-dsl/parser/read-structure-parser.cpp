#include "adapter-dsl/parser/read-structure-parser.h"

#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include <xoos/error/error.h>

namespace xoos::demux::dsl {

const ParserRegistry& ParserRegistry::Instance() {
  // Function-local Meyers singleton, not a mutable global; inline variables do not apply.
  static const ParserRegistry instance = [] {  // NOSONAR(S6018)
    ParserRegistry registry;
    RegisterBuiltinParsers(registry);
    return registry;
  }();
  return instance;
}

void ParserRegistry::Register(std::shared_ptr<ReadStructureParser> parser) {
  if (parser == nullptr) {
    throw error::Error("Cannot register a null read-structure parser");
  }
  const std::string_view syntax = parser->Syntax();
  for (auto& [name, existing] : _parsers) {
    if (name == syntax) {
      existing = std::move(parser);
      return;
    }
  }
  _parsers.emplace_back(std::string{syntax}, std::move(parser));
}

const ReadStructureParser& ParserRegistry::Get(std::string_view syntax) const {
  for (const auto& [name, parser] : _parsers) {
    if (name == syntax) {
      return *parser;
    }
  }
  std::string known;
  for (const auto& [name, parser] : _parsers) {
    (void)parser;
    if (!known.empty()) {
      known += ", ";
    }
    known += name;
  }
  throw error::Error("Unknown read-structure syntax '{}'. Known syntaxes: {}", syntax, known);
}

CompiledAdapter ParseReadStructure(const std::string_view syntax, const std::string_view text) {
  return ParserRegistry::Instance().Get(syntax).Parse(text);
}

}  // namespace xoos::demux::dsl
