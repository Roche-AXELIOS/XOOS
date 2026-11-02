#include <memory>
#include <string_view>

#include "adapter-dsl/dsl/dsl-compiler.h"
#include "adapter-dsl/model/adapter-model.h"
#include "adapter-dsl/parser/read-structure-parser.h"

namespace xoos::demux::dsl {

namespace {

/**
 * The shipping DSL parser, adapting the free-function @ref CompileDsl to the parser interface.
 * This is the single place the "dsl" syntax name binds to the DSL implementation; swapping the
 * DSL frontend (or adding a second syntax) is a local change here, not at any call site.
 */
class DslParser final : public ReadStructureParser {
 public:
  std::string_view Syntax() const override {
    return kDefaultSyntax;
  }

  CompiledAdapter Parse(const std::string_view text) const override {
    return CompileDsl(text);
  }
};

}  // namespace

/**
 * Register every parser shipped with demux. Called once, lazily, by ParserRegistry::Instance().
 * Add a new authoring syntax by implementing ReadStructureParser and registering it here.
 */
void RegisterBuiltinParsers(ParserRegistry& registry) {
  registry.Register(std::make_shared<DslParser>());
}

}  // namespace xoos::demux::dsl
