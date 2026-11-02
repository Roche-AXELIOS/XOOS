#pragma once

#include <cstddef>
#include <string>

#include <xoos/cli/cli.h>
#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "core/read-collapser-options.h"

namespace xoos::read_collapser {

constexpr u32 kU8Limit = 255;

constexpr s32 kMaxReadsPerCluster = 50;
constexpr s32 kMinReadsPerCluster = 2;
constexpr s32 kMinCompressionLevel = 1;
constexpr s32 kMaxCompressionLevel = 9;

constexpr s32 kDefaultPadding = 200;
constexpr s32 kDefaultWiggleRoom = 0;
constexpr s32 kDefaultWiggleRoomPartial = 0;
constexpr s32 kDefaultBatchSize = 300'000;
constexpr u32 kDefaultRegionSize = 10'000;
const u16 kDefaultExcludeFlag = BAM_FSUPPLEMENTARY | BAM_FSECONDARY;
const u16 kDefaultExcludeFlagConsensus = BAM_FSUPPLEMENTARY | BAM_FSECONDARY;
const u16 kDefaultExcludeFlagMarkdup = BAM_FSECONDARY;
constexpr f64 kDefaultConsensusThreshold = 0.5;
constexpr f64 kDefaultConsensusGapThreshold = 0.5;

const std::string kOptGroupNamePresetOptions = "Preset Options";
const std::string kOptGroupNameInputOptions = "Input Options";
const std::string kOptGroupNameOutputOptions = "Output Options";
const std::string kOptGroupNameRegionOptions = "Region Options";
const std::string kOptGroupNameReadFilterOptions = "Read Filtering Options";
const std::string kOptGroupNameClusterOptions = "Cluster Options";
const std::string kOptGroupNameConsensusOptions = "Consensus Options";
const std::string kOptGroupNamePerformanceOptions = "Performance Options";
const std::string kOptGroupNameHiddenOptions = "";  // NOLINT

struct PresetHash {
  size_t operator()(const ReadCollapserPresets& presets) const;
};

using PresetDefaultsMap = std::unordered_map<std::string, std::string>;
using PresetsMap = std::unordered_map<ReadCollapserPresets, PresetDefaultsMap, PresetHash>;

// markdup presets
const PresetsMap kMarkdupPresetsMap = {
    {ReadCollapserPresets::kWgsDuplex, {{"--cluster-by-strand", "true"}, {"--wiggle-room", "0"}}},
    {ReadCollapserPresets::kWgsSimplex,
     {{"--make-clusters-of-partial-reads-only", "true"}, {"--wiggle-room", "2"}, {"--wiggle-room-partial", "0"}}},
    {ReadCollapserPresets::kRnaBulk, {{"--cluster-by-strand", "true"}, {"--wiggle-room", "0"}}},
};

// consensus presets
const PresetsMap kConsensusPresetsMap = {
    {ReadCollapserPresets::kWgsDuplex, {{"--cluster-by-strand", "true"}, {"--wiggle-room", "0"}}},
    {ReadCollapserPresets::kWgsDuplexCfdna,
     {{"--cluster-by-strand", "true"},
      {"--wiggle-room", "2"},
      {"--duplex-library-type", "parent-parent"},
      {"--min-cluster-size", "1"},
      {"--min-trim-read-support", "1"},
      {"--max-discordant-duplex-error-percentage", "5"}}},
    {ReadCollapserPresets::kTeDuplex,
     {{"--cluster-by-umi", "true"},
      {"--wiggle-room", "2"},
      {"--wiggle-room-partial", "0"},
      {"--duplex-library-type", "parent-daughter"},
      {"--min-same-strand-cluster-size", "3"}}},
    {ReadCollapserPresets::kWgsSimplex,
     {{"--make-clusters-of-partial-reads-only", "true"},
      {"--wiggle-room", "2"},
      {"--wiggle-room-partial", "0"},
      {"--min-cluster-size", "1"},
      {"--min-trim-read-support", "1"}}},
    {ReadCollapserPresets::kTeSimplex,
     {{"--cluster-by-umi", "true"}, {"--wiggle-room", "2"}, {"--wiggle-room-partial", "0"}}},
};

void AddPresetOption(CLI::App* app, const PresetsMap& preset_map, const std::string& group_name);

void AddCommonInputOptions(CLI::App* app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddCommonOutputOptions(CLI::App* app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddCommonRegionOptions(CLI::App* app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddCommonReadFilteringOptions(CLI::App* app,
                                   const ReadCollapserOptionsPtr& options,
                                   const std::string& group_name);

void AddMarkdupReadFilterOptions(CLI::App* app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddConsensusReadFilterOptions(CLI::App* app,
                                   const ReadCollapserOptionsPtr& options,
                                   const std::string& group_name);

void AddCommonClusterOptions(CLI::App* app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddCommonPerformanceOptions(CLI::App* app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddConsensusOptions(CLI::App* app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddMarkdupOutputOptions(CLI::App* app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddConsensusOutputOptions(CLI::App* app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddConsensusClusterOptions(CLI::App* app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddCommonHiddenOptions(CLI::App* app, const ReadCollapserOptionsPtr& options);

void ValidateConsensusOptions(const CLI::App* app, const ReadCollapserOptionsPtr& options);

void ValidateMarkdupOptions(const CLI::App* app, const ReadCollapserOptionsPtr& options);

}  // namespace xoos::read_collapser
