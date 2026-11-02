#include "adapter-dsl/authored/authored-json.h"

#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include <nlohmann/json.hpp>

#include <xoos/error/error.h>
#include <xoos/vfs/vfs-iterator.h>

#include "adapter-dsl/model/adapter-entry-json.h"
#include "adapter-dsl/model/adapter-model-constants.h"
#include "adapter-dsl/model/adapter-model.h"
#include "adapter-dsl/parser/read-structure-parser.h"

namespace xoos::demux::dsl {

using json = nlohmann::json;  // NOLINT

namespace {

/**
 * Collect fixed literal sequences from the DSL-compiled read structure. A `type: fixed` side-car
 * entry must be keyed by one of these literals.
 */
std::set<std::string, std::less<>> CollectFixedSequences(const CompiledAdapter& adapter) {
  std::set<std::string, std::less<>> sequences;
  for (const Token& token : adapter.read_structure) {
    std::visit(
        [&sequences]<typename T>([[maybe_unused]] const T& value) {
          if constexpr (std::is_same_v<T, FixedToken>) {
            sequences.insert(value.sequence);
          } else if constexpr (std::is_same_v<T, ArmToken>) {
            for (const LeafToken& leaf : value.tokens) {
              if (const auto* const fixed = std::get_if<FixedToken>(&leaf)) {
                sequences.insert(fixed->sequence);
              }
            }
          }
        },
        token);
  }
  return sequences;
}

/**
 * Overlay the authored `adapters` block onto the bare skeleton produced by the DSL compiler.
 *
 * The DSL fixes the set of ids and their sid/umi types. The authored block supplies the data the
 * DSL cannot express: the FASTA `sequence` and the v1 legacy params for each id. Extra entries are
 * allowed only for fixed literal side-cars keyed by an uppercase fixed sequence in the read
 * structure. We require both sides to agree on the id set and declared type, so a typo in either
 * surfaces loudly.
 */
void OverlayAdapters(const std::string& name, CompiledAdapter& adapter, const json& authored_adapters) {
  for (auto& [base, entry] : adapter.adapters) {
    if (!authored_adapters.contains(base)) {
      throw error::Error(
          "Adapter '{}': id '{}' from the DSL has no matching entry in the authored adapters block", name, base);
    }
    // Keep the DSL-derived type; overlay sequence + legacy params (and verify any declared type).
    OverlayEntryFromJson(base, entry, authored_adapters.at(base));
  }
  const std::set<std::string, std::less<>> fixed_sequences = CollectFixedSequences(adapter);
  for (const auto& item : authored_adapters.items()) {
    const std::string& base = item.key();
    if (adapter.Find(base) != nullptr) {
      // already overlaid above as a DSL id
      continue;
    }
    AdapterEntry entry = EntryFromJson(base, item.value());
    if (entry.type != EntryType::kFixed) {
      throw error::Error("Adapter '{}': adapters block has entry '{}' that the DSL never references", name, base);
    }
    if (!fixed_sequences.contains(base)) {
      throw error::Error(
          "Adapter '{}': fixed adapters entry '{}' does not match any fixed literal in the read structure", name, base);
    }
    adapter.adapters.emplace_back(base, std::move(entry));
  }
}

CompiledAdapter CompileAuthoredAdapter(const std::string& name, const json& j) {
  if (!j.contains(kDslKey)) {
    throw error::Error("Authored adapter '{}' must contain a '{}' string", name, kDslKey);
  }
  // The authoring syntax is selectable via an optional `syntax` field; absent, it defaults to the
  // DSL. This is the only place a syntax name meets an authored string (see ParserRegistry).
  const std::string syntax = j.value(std::string{kSyntaxKey}, std::string{kDefaultSyntax});
  CompiledAdapter adapter = ParseReadStructure(syntax, j.at(kDslKey).get<std::string>());
  if (!j.contains(kTypeKey)) {
    throw error::Error("Authored adapter '{}' must declare a '{}' (the design family)", name, kTypeKey);
  }
  adapter.type = ParseDesignType(j.at(std::string{kTypeKey}).get<std::string>());
  const json authored_adapters = j.contains(kAdaptersKey) ? j.at(kAdaptersKey) : json::object();
  OverlayAdapters(name, adapter, authored_adapters);
  return adapter;
}

}  // namespace

AdapterBundle FromAuthoredJson(const std::string& json_text) {
  json parsed;
  try {
    parsed = json::parse(json_text, /*callback=*/nullptr, /*allow_exceptions=*/true, /*ignore_comments=*/true);
  } catch (const json::exception& ex) {
    throw error::Error("Authored adapter JSON cannot be parsed: {}", ex.what());
  }
  if (!parsed.is_object()) {
    throw error::Error("Authored adapter bundle must be a JSON object");
  }
  // Two accepted shapes:
  //   1. Wrapper (preferred): {"adapter_designs": { "<name>": {<adapter>}, ... }}
  //   2. Flat map: {"<name>": {<adapter>}, ...}.
  AdapterBundle bundle;
  if (parsed.contains(kAdapterDesignsKey)) {
    const json& designs = parsed.at(kAdapterDesignsKey);
    if (!designs.is_object()) {
      throw error::Error("'{}' must be a JSON object mapping names to adapters", kAdapterDesignsKey);
    }
    for (const auto& [name, value] : designs.items()) {
      bundle.adapters.emplace_back(name, CompileAuthoredAdapter(name, value));
    }
  } else {
    for (const auto& [name, value] : parsed.items()) {
      bundle.adapters.emplace_back(name, CompileAuthoredAdapter(name, value));
    }
  }
  ValidateBundle(bundle);
  return bundle;
}

AdapterBundle FromAuthoredJsonFile(const fs::path& path) {
  std::ifstream stream{path};
  if (!stream.is_open()) {
    throw error::Error("Authored adapter JSON file '{}' could not be opened", path.string());
  }
  std::ostringstream buffer;
  buffer << stream.rdbuf();
  return FromAuthoredJson(buffer.str());
}

AdapterBundle FromAuthoredJsonBundle(const vfs::VirtualFilesystemPtr& vfs, const fs::path& path) {
  const vfs::VirtualFileHandlePtr fh = vfs->Open(path);
  if (fh == nullptr) {
    throw error::Error("Adapter config '{}' does not exist in bundle '{}'", path.string(), vfs->GetName());
  }
  const std::string text{vfs::begin(fh), vfs::end(fh)};
  return FromAuthoredJson(text);
}

}  // namespace xoos::demux::dsl
