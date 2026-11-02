#pragma once

#include <string>

#include <nlohmann/json.hpp>

#include "adapter-dsl/model/adapter-model.h"

namespace xoos::demux::dsl {

/// @brief Serialize a Wiggle window {left, right} to JSON. Shared with the fixed-token codec.
nlohmann::json WiggleToJson(const Wiggle& w);

/// @brief Parse a Wiggle window {left, right} from JSON. Shared with the fixed-token codec.
Wiggle WiggleFromJson(const nlohmann::json& j);

/**
 * @brief Serialize an AdapterEntry to its JSON object form (type, sequence, v1 legacy params).
 * Shared by the compiled and authored JSON readers/writers so the `adapters` shape stays identical
 * in both representations.
 */
nlohmann::json EntryToJson(const AdapterEntry& entry);

/**
 * @brief Parse an AdapterEntry from its JSON object form. @p base is used only for error context.
 * @throws xoos::Error if the required `type` field is missing or invalid.
 */
AdapterEntry EntryFromJson(const std::string& base, const nlohmann::json& j);

/**
 * @brief Overlay the JSON fields onto @p entry: set sequence and the v1 legacy params.
 * The entry's existing `type` (from the DSL) is retained. If the JSON carries a `type` it must
 * match @p entry.type, otherwise an xoos::Error is thrown. Used by the authored-form overlay
 * where type is inferred from the DSL and the adapters block supplies the rest.
 */
void OverlayEntryFromJson(const std::string& base, AdapterEntry& entry, const nlohmann::json& j);

}  // namespace xoos::demux::dsl
