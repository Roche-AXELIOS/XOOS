#include "adapter-dsl/model/compiled-json.h"

#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <variant>

#include <nlohmann/json.hpp>

#include <xoos/error/error.h>
#include <xoos/types/int.h>

#include "adapter-dsl/model/adapter-entry-json.h"
#include "adapter-dsl/model/adapter-model-constants.h"
#include "adapter-dsl/model/adapter-model.h"

namespace xoos::demux::dsl {

using json = nlohmann::json;  // NOLINT

namespace {

/// Default indentation width (in spaces) for serialized compiled-form JSON.
constexpr s32 kDefaultJsonIndent = 2;

// --- read_structure <-> JSON ---

// Per the design spec, every leaf token in `read_structure` serializes to a bare string:
// a fixed leaf is its UPPERCASE sequence, an id leaf is `base` + orientation suffix, and the insert
// keyword is `insert`(+/-). Optional `type: fixed` side-cars live in `adapters`, not on the token.
// A duplex arm serializes to a nested array.
json LeafToJson(const LeafToken& leaf) {
  return LeafToCompiledString(leaf);
}

LeafToken LeafFromJson(const std::string& name, const json& j) {
  if (j.is_string()) {
    return LeafFromCompiledString(j.get<std::string>());
  }
  throw error::Error("Adapter '{}': a leaf read_structure token must be a string", name);
}

json TokenToJson(const Token& token) {
  if (const auto* const arm = std::get_if<ArmToken>(&token)) {
    json array = json::array();
    for (const LeafToken& leaf : arm->tokens) {
      array.push_back(LeafToJson(leaf));
    }
    return array;
  }
  // Leaf token: fixed / id / insert.
  if (const auto* const fixed = std::get_if<FixedToken>(&token)) {
    return LeafToJson(LeafToken{*fixed});
  }
  if (const auto* const id = std::get_if<IdToken>(&token)) {
    return LeafToJson(LeafToken{*id});
  }
  return LeafToJson(LeafToken{std::get<InsertToken>(token)});
}

Token TokenFromJson(const std::string& name, const json& j) {
  if (j.is_string()) {
    // NOLINTNEXTLINE(whitespace/newline)
    return std::visit([](auto leaf) { return Token{leaf}; }, LeafFromJson(name, j));
  }
  if (j.is_array()) {
    ArmToken arm;
    for (const auto& element : j) {
      if (!element.is_string()) {
        throw error::Error("Adapter '{}': a duplex arm may only contain string tokens", name);
      }
      arm.tokens.push_back(LeafFromJson(name, element));
    }
    return arm;
  }
  throw error::Error("Adapter '{}': read_structure token must be a string or an array", name);
}

// --- CompiledAdapter <-> JSON ---

json AdapterToJson(const CompiledAdapter& adapter) {
  json j = json::object();
  j[std::string{kTypeKey}] = ToString(adapter.type);
  json structure = json::array();
  for (const Token& token : adapter.read_structure) {
    structure.push_back(TokenToJson(token));
  }
  j[std::string{kReadStructureKey}] = std::move(structure);
  json entries = json::object();
  for (const auto& [base, entry] : adapter.adapters) {
    entries[base] = EntryToJson(entry);
  }
  j[std::string{kAdaptersKey}] = std::move(entries);
  return j;
}

CompiledAdapter AdapterFromJson(const std::string& name, const json& j) {
  CompiledAdapter adapter;
  if (!j.contains(kTypeKey)) {
    throw error::Error("Adapter '{}' must declare a '{}' (the design family)", name, kTypeKey);
  }
  adapter.type = ParseDesignType(j.at(std::string{kTypeKey}).get<std::string>());
  if (!j.contains(kReadStructureKey)) {
    throw error::Error("Adapter '{}' must contain a '{}' array", name, kReadStructureKey);
  }
  for (const auto& token : j.at(kReadStructureKey)) {
    adapter.read_structure.push_back(TokenFromJson(name, token));
  }
  if (j.contains(kAdaptersKey)) {
    // Preserve entry order as written in the JSON object.
    const auto& adapter_entries = j.at(kAdaptersKey);
    for (const auto& [base, entry] : adapter_entries.items()) {
      adapter.adapters.emplace_back(base, EntryFromJson(base, entry));
    }
  }
  return adapter;
}

}  // namespace

std::string CompiledToJson(const AdapterBundle& bundle, const s32 indent) {
  json j = json::object();
  for (const auto& [name, adapter] : bundle.adapters) {
    j[name] = AdapterToJson(adapter);
  }
  return j.dump(indent);
}

std::string CompiledToJson(const AdapterBundle& bundle) {
  return CompiledToJson(bundle, kDefaultJsonIndent);
}

AdapterBundle CompiledFromJson(const std::string& json_text) {
  json parsed;
  try {
    parsed = json::parse(json_text);
  } catch (const json::exception& ex) {
    throw error::Error("Compiled adapter JSON cannot be parsed: {}", ex.what());
  }
  if (!parsed.is_object()) {
    throw error::Error("Compiled adapter bundle must be a JSON object mapping names to adapters");
  }
  AdapterBundle bundle;
  for (const auto& [name, value] : parsed.items()) {
    bundle.adapters.emplace_back(name, AdapterFromJson(name, value));
  }
  ValidateBundle(bundle);
  return bundle;
}

AdapterBundle CompiledFromJsonFile(const fs::path& path) {
  std::ifstream stream{path};
  if (!stream.is_open()) {
    throw error::Error("Compiled adapter JSON file '{}' could not be opened", path.string());
  }
  std::ostringstream buffer;
  buffer << stream.rdbuf();
  return CompiledFromJson(buffer.str());
}

}  // namespace xoos::demux::dsl
