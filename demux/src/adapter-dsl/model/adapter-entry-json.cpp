#include "adapter-dsl/model/adapter-entry-json.h"

#include <algorithm>
#include <initializer_list>
#include <optional>
#include <string>
#include <string_view>

#include <nlohmann/json.hpp>

#include <xoos/error/error.h>
#include <xoos/types/int.h>

#include "adapter-dsl/model/adapter-model-constants.h"
#include "adapter-dsl/model/adapter-model.h"

namespace xoos::demux::dsl {

using json = nlohmann::json;  // NOLINT

json WiggleToJson(const Wiggle& w) {
  return json{{std::string{kWiggleLeftKey}, w.left}, {std::string{kWiggleRightKey}, w.right}};
}

Wiggle WiggleFromJson(const nlohmann::json& j) {
  Wiggle w;
  w.left = j.value(std::string{kWiggleLeftKey}, s64{0});
  w.right = j.value(std::string{kWiggleRightKey}, s64{0});
  return w;
}

namespace {

bool ContainsAny(const json& j, const std::initializer_list<std::string_view> keys) {
  return std::ranges::any_of(keys, [&j](const std::string_view key) { return j.contains(key); });
}

void SetOptionalInt(json& j, const std::string_view key, const std::optional<s64>& value) {
  if (value.has_value()) {
    j[std::string{key}] = value.value();
  }
}

std::optional<s64> GetOptionalInt(const json& j, std::string_view key) {
  return j.contains(key) ? std::make_optional(j.at(key).get<s64>()) : std::nullopt;
}

/// Read sequence + v1 legacy params from @p j into @p entry (leaves entry.type untouched).
void ReadSequenceAndParams(AdapterEntry& entry, const json& j) {
  entry.sequence = j.value(std::string{kSequenceKey}, std::string{});
  entry.max_edit_distance = GetOptionalInt(j, kMaxEditDistanceKey);
  if (j.contains(kWiggleFwKey)) {
    entry.wiggle_forward = WiggleFromJson(j.at(kWiggleFwKey));
  }
  if (j.contains(kWiggleRcKey)) {
    entry.wiggle_reverse = WiggleFromJson(j.at(kWiggleRcKey));
  }
  entry.bait_length_forward = GetOptionalInt(j, kBaitLengthFwKey);
  entry.bait_length_reverse = GetOptionalInt(j, kBaitLengthRcKey);
  entry.additional_prefix_length = GetOptionalInt(j, kAdditionalPrefixLengthKey);
  entry.additional_suffix_length = GetOptionalInt(j, kAdditionalSuffixLengthKey);
}

/// Validate and normalize a fixed side-car. Fixed entries are metadata only; no matcher knobs.
void ReadFixedLiteral(AdapterEntry& entry, const std::string& base, const json& j) {
  if (!IsUpperAcgt(base)) {
    throw error::Error("fixed adapters entry '{}' must be keyed by an uppercase ACGT literal", base);
  }
  if (ContainsAny(j,
                  {kMaxEditDistanceKey,
                   kWiggleFwKey,
                   kWiggleRcKey,
                   kBaitLengthFwKey,
                   kBaitLengthRcKey,
                   kAdditionalPrefixLengthKey,
                   kAdditionalSuffixLengthKey})) {
    throw error::Error("fixed adapters entry '{}' must not contain matcher tuning fields", base);
  }
  const std::string sequence = j.value(std::string{kSequenceKey}, base);
  if (sequence != base) {
    throw error::Error("fixed adapters entry '{}' sequence '{}' must match its key", base, sequence);
  }
  entry.sequence = base;
}

}  // namespace

json EntryToJson(const AdapterEntry& entry) {
  json j = json::object();
  j[std::string{kEntryTypeKey}] = ToString(entry.type);
  j[std::string{kSequenceKey}] = entry.sequence;
  SetOptionalInt(j, kMaxEditDistanceKey, entry.max_edit_distance);
  if (entry.wiggle_forward.has_value()) {
    j[std::string{kWiggleFwKey}] = WiggleToJson(entry.wiggle_forward.value());
  }
  if (entry.wiggle_reverse.has_value()) {
    j[std::string{kWiggleRcKey}] = WiggleToJson(entry.wiggle_reverse.value());
  }
  SetOptionalInt(j, kBaitLengthFwKey, entry.bait_length_forward);
  SetOptionalInt(j, kBaitLengthRcKey, entry.bait_length_reverse);
  SetOptionalInt(j, kAdditionalPrefixLengthKey, entry.additional_prefix_length);
  SetOptionalInt(j, kAdditionalSuffixLengthKey, entry.additional_suffix_length);
  return j;
}

AdapterEntry EntryFromJson(const std::string& base, const nlohmann::json& j) {
  AdapterEntry entry;
  if (!j.contains(kEntryTypeKey)) {
    throw error::Error("adapters entry '{}' must declare a 'type' (sid / umi / fixed)", base);
  }
  entry.type = ParseEntryType(j.at(kEntryTypeKey).get<std::string>());
  if (entry.type == EntryType::kFixed) {
    ReadFixedLiteral(entry, base, j);
  } else {
    ReadSequenceAndParams(entry, j);
  }
  return entry;
}

void OverlayEntryFromJson(const std::string& base, AdapterEntry& entry, const nlohmann::json& j) {
  // type is optional in the authored form (inferred from the DSL); if present, it must agree.
  if (j.contains(kEntryTypeKey)) {
    const EntryType declared = ParseEntryType(j.at(kEntryTypeKey).get<std::string>());
    if (declared != entry.type) {
      throw error::Error("adapters entry '{}' is declared as '{}' in the DSL but '{}' in the adapters block",
                         base,
                         ToString(entry.type),
                         ToString(declared));
    }
  }
  ReadSequenceAndParams(entry, j);
}

}  // namespace xoos::demux::dsl
