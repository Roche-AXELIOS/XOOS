#include "adapter-design/adapter-design.h"

#include <array>
#include <format>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

#include <nlohmann/json.hpp>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/vfs/vfs.h>

#include "adapter-design/adapter-design-constants.h"
#include "adapter-dsl/assemble/adapter-design-assembler.h"
#include "adapter-dsl/authored/authored-json.h"
#include "adapter-dsl/model/adapter-model.h"
#include "sequence/sequence-enums.h"

using json = nlohmann::json;  // NOLINT

namespace {

/// Length of the trailing separator ", " appended after each enum value name.
constexpr std::size_t kSeparatorLength = 2;

/**
 * Lookup a JSON string value in a mapping table and return the corresponding enum.
 * Throws if the value is not found, listing all valid options.
 */
template <typename EnumType, std::size_t N>
EnumType StrictEnumLookup(const json& j,
                          const std::array<std::pair<EnumType, std::string_view>, N>& mapping,
                          std::string_view type_name) {
  const auto str = j.get<std::string>();
  for (const auto& [value, name] : mapping) {
    if (name == str) {
      return value;
    }
  }
  std::string valid;
  for (const auto& [_, name] : mapping) {
    valid += std::format("'{}', ", name);
  }
  if (valid.size() >= kSeparatorLength) {
    valid.resize(valid.size() - kSeparatorLength);
  }
  throw xoos::error::Error("Unknown {} '{}', valid values: {}", type_name, str, valid);
}

/**
 * Lookup an enum value in a mapping table and return the corresponding JSON string.
 * Throws if the value is not found.
 */
template <typename EnumType, std::size_t N>
std::string_view StrictEnumToString(EnumType e,
                                    const std::array<std::pair<EnumType, std::string_view>, N>& mapping,
                                    std::string_view type_name) {
  for (const auto& [value, name] : mapping) {
    if (value == e) {
      return name;
    }
  }
  throw xoos::error::Error("Unknown {} enum value {}", type_name, static_cast<xoos::s32>(e));
}

}  // namespace

namespace xoos::demux {

// --- Enum mapping tables ---

static constexpr auto kMatcherKindMapping = std::to_array<std::pair<MatcherKind, std::string_view>>({
    {MatcherKind::kUmi, "umi"},
    {MatcherKind::kSid, "sid"},
});

constexpr std::string_view kFixedComponentType = "fixed";
constexpr std::string_view kMaxWiggleLeftKey = "max_wiggle_left";
constexpr std::string_view kMaxWiggleRightKey = "max_wiggle_right";

static constexpr auto kSeqTransformMapping = std::to_array<std::pair<SeqTransform, std::string_view>>({
    {SeqTransform::kNone, "none"},
    {SeqTransform::kReverse, "reverse"},
    {SeqTransform::kReverseComplement, "reverse_complement"},
});

static constexpr auto kAdapterTypeMapping = std::to_array<std::pair<AdapterType, std::string_view>>({
    {AdapterType::kDuplex, "duplex"},
    {AdapterType::kDuplexUMI, "duplex_umi"},
    {AdapterType::kDuplexStem, "duplex_stem"},
    {AdapterType::kYsu, "YSU"},
    {AdapterType::kYs, "YS"},
    {AdapterType::kSimplex, "simplex"},
    {AdapterType::kSimplex10x, "simplex-10x"},
});

// --- Enum to_json / from_json ---

void to_json(json& j, const MatcherKind& e) {  // NOSONAR(S100)
  j = StrictEnumToString(e, kMatcherKindMapping, "MatcherKind");
}

void from_json(const json& j, MatcherKind& e) {  // NOSONAR(S100)
  e = StrictEnumLookup(j, kMatcherKindMapping, "MatcherKind");
}

void to_json(json& j, const SeqTransform& e) {  // NOSONAR(S100)
  j = StrictEnumToString(e, kSeqTransformMapping, "SeqTransform");
}

void from_json(const json& j, SeqTransform& e) {  // NOSONAR(S100)
  e = StrictEnumLookup(j, kSeqTransformMapping, "SeqTransform");
}

void to_json(json& j, const AdapterType& e) {  // NOSONAR(S100)
  j = StrictEnumToString(e, kAdapterTypeMapping, "AdapterType");
}

void from_json(const json& j, AdapterType& e) {  // NOSONAR(S100)
  e = StrictEnumLookup(j, kAdapterTypeMapping, "AdapterType");
}

// --- Struct to_json / from_json ---

void to_json(json& j, const SearchWindow& d) {  // NOSONAR(S100)
  j = json{{kMaxWiggleLeftKey, d.max_wiggle_left}, {kMaxWiggleRightKey, d.max_wiggle_right}};
}

void from_json(const json& j, SearchWindow& d) {  // NOSONAR(S100)
  j.at(kMaxWiggleLeftKey).get_to(d.max_wiggle_left);
  j.at(kMaxWiggleRightKey).get_to(d.max_wiggle_right);
}

void to_json(json& j, const FixedComponentDefinition& d) {  // NOSONAR(S100)
  j = json{{kTypeKey, kFixedComponentType}, {kSequenceKey, d.sequence}};
}

void from_json(const json& j, FixedComponentDefinition& d) {  // NOSONAR(S100)
  const auto type = j.at(kTypeKey).get<std::string>();
  if (type != kFixedComponentType) {
    throw error::Error("Expected fixed adapter component, found type '{}'", type);
  }
  if (j.contains(kSequencesKey)) {
    throw error::Error("Fixed adapter component must not contain matcher pool '{}'", kSequencesKey);
  }
  j.at(kSequenceKey).get_to(d.sequence);
  if (d.sequence.empty()) {
    throw error::Error("Fixed adapter component '{}' must not be empty", kSequenceKey);
  }
}

void to_json(json& j, const MatcherComponentDefinition& d) {  // NOSONAR(S100)
  j = json{{kTypeKey, d.kind},
           {kSequencesKey, d.sequences},
           {kMaxEditDistanceKey, d.max_edit_distance},
           {kTransformKey, d.transform},
           {kMaxWiggleLeftKey, d.search.max_wiggle_left},
           {kMaxWiggleRightKey, d.search.max_wiggle_right}};
  if (!d.prefix.empty()) {
    j[kPrefixKey] = d.prefix;
  }
  if (!d.suffix.empty()) {
    j[kSuffixKey] = d.suffix;
  }
  if (!d.bait.empty()) {
    j[kBaitKey] = d.bait;
  }
}

void from_json(const json& j, MatcherComponentDefinition& d) {  // NOSONAR(S100)
  j.at(kTypeKey).get_to(d.kind);
  j.at(kSequencesKey).get_to(d.sequences);
  if (d.sequences.empty() && d.kind != MatcherKind::kSid) {
    throw error::Error("Matcher adapter component '{}' must not be empty for non-SID matchers", kSequencesKey);
  }
  j.at(kMaxEditDistanceKey).get_to(d.max_edit_distance);
  j.at(kTransformKey).get_to(d.transform);
  d.prefix = j.value(kPrefixKey, std::string{});
  d.suffix = j.value(kSuffixKey, std::string{});
  d.search.max_wiggle_left = j.at(kMaxWiggleLeftKey).get<s32>();
  d.search.max_wiggle_right = j.at(kMaxWiggleRightKey).get<s32>();
  d.bait = j.value(kBaitKey, std::string{});
}

void to_json(json& j, const AdapterComponentDefinition& d) {  // NOSONAR(S100)
  std::visit([&j](const auto& component) { to_json(j, component); }, d);
}

void from_json(const json& j, AdapterComponentDefinition& d) {  // NOSONAR(S100)
  const auto type = j.at(kTypeKey).get<std::string>();
  if (type == kFixedComponentType) {
    d = j.get<FixedComponentDefinition>();
  } else {
    d = j.get<MatcherComponentDefinition>();
  }
}

void to_json(json& j, const AdapterDesign& d) {  // NOSONAR(S100)
  j = json{{kNameKey, d.name}, {kTypeKey, d.type}, {kComponentsKey, d.components}};
}

void from_json(const json& j, AdapterDesign& d) {  // NOSONAR(S100)
  j.at(kNameKey).get_to(d.name);
  j.at(kTypeKey).get_to(d.type);
  j.at(kComponentsKey).get_to(d.components);
}

std::string FormatJson(const AdapterComponentDefinition& component) {
  json j;
  to_json(j, component);
  return j.dump();
}

AdapterDesign LoadAdapterDesign(const fs::path& adapter_design_bundle, const std::string& adapter_design_name) {
  vfs::VirtualFilesystemPtr vfs = vfs::Open(adapter_design_bundle);
  const dsl::AdapterBundle bundle = dsl::FromAuthoredJsonBundle(vfs, fs::path{kAdapterConfigFilename});
  Logging::Info("Loading adapter design '{}' from bundle '{}'", adapter_design_name, adapter_design_bundle.string());
  const dsl::CompiledAdapter& adapter = bundle.At(adapter_design_name);
  return dsl::ToAdapterDesign(adapter_design_name, adapter);
}
}  // namespace xoos::demux
