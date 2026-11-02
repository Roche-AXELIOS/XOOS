#pragma once
#include <string_view>

namespace xoos::demux {

constexpr std::string_view kNameKey = "name";
constexpr std::string_view kTypeKey = "type";
constexpr std::string_view kComponentsKey = "components";
constexpr std::string_view kBaitKey = "bait";

constexpr std::string_view kSequenceKey = "sequence";
constexpr std::string_view kSequencesKey = "sequences";
constexpr std::string_view kMaxEditDistanceKey = "max_edit_distance";
constexpr std::string_view kTransformKey = "transform";
constexpr std::string_view kPrefixKey = "prefix";
constexpr std::string_view kSuffixKey = "suffix";

/// @brief Name of the DSL adapter config file inside the adapter design bundle.
constexpr std::string_view kAdapterConfigFilename = "adapter-config.json5";

}  // namespace xoos::demux
