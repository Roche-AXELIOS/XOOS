#pragma once

#include <xoos/cli/cli.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "core/read-collapser-options.h"

namespace xoos::read_collapser {

constexpr u32 kU8Limit = 255;

constexpr s32 kMaxReadsPerCluster = 50;
constexpr s32 kMinReadsPerCluster = 2;
constexpr s32 kMinCompressionLevel = 1;
constexpr s32 kMaxCompressionLevel = 9;

constexpr s32 kDefaultPadding = 200;
constexpr s32 kDefaultWiggleRoom = 4;
constexpr s32 kDefaultWiggleRoomPartial = 2;
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
const std::string kOptGroupNameReadFilterOptions = "Read-level Filter Options";
const std::string kOptGroupNameClusterOptions = "Cluster Options";
const std::string kOptGroupNameConsensusOptions = "Consensus Options";
const std::string kOptGroupNameConsensusDebugOptions = "Consensus Debug Options";
const std::string kOptGroupNamePerformanceOptions = "Performance Options";

struct PresetHash {
  size_t operator()(const ReadCollapserPresets& presets) const;
};

using PresetDefaultsMap = std::unordered_map<std::string, std::string>;
using PresetsMap = std::unordered_map<ReadCollapserPresets, PresetDefaultsMap, PresetHash>;

// markdup presets
const PresetsMap kMarkdupPresetsMap = {
    {ReadCollapserPresets::kWgsDuplex, {{"--cluster-by-strand", "true"}, {"--wiggle-room", "0"}}},
    {ReadCollapserPresets::kWgsSimplex,
     {{"--cluster-by-strand", "true"},
      {"--cluster-by-umi", "true"},
      {"--make-clusters-of-partial-reads-only", "true"},
      {"--wiggle-room", "0"}}},
    {ReadCollapserPresets::kRnaBulk, {{"--cluster-by-strand", "true"}, {"--wiggle-room", "0"}}},
};

// consensus presets
const PresetsMap kConsensusPresetsMap = {
    {ReadCollapserPresets::kWgsDuplex, {{"--cluster-by-strand", "true"}, {"--wiggle-room", "0"}}},
    {ReadCollapserPresets::kWgsDuplexMrd,
     {{"--cluster-by-strand", "true"},
      {"--wiggle-room", "2"},
      {"--duplex-library-type", "parent-parent"},
      {"--min-cluster-size", "1"},
      {"--min-trim-read-support", "1"},
      {"--max-discordant-duplex-error-percentage", "5"}}},
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
    {ReadCollapserPresets::kTeSimplex,
     {{"--cluster-by-umi", "true"}, {"--wiggle-room", "2"}, {"--wiggle-room-partial", "0"}}},
};

void AddPresetOption(cli::AppPtr app, const PresetsMap& preset_map, const std::string& group_name);

void AddCommonInputOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddCommonOutputOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddCommonReadFilteringOptions(cli::AppPtr app,
                                   const ReadCollapserOptionsPtr& options,
                                   const std::string& group_name);

void AddMarkdupReadFilterOptions(cli::AppPtr app,
                                 const ReadCollapserOptionsPtr& options,
                                 const std::string& group_name);

void AddConsensusReadFilterOptions(cli::AppPtr app,
                                   const ReadCollapserOptionsPtr& options,
                                   const std::string& group_name);

void AddCommonClusterOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddCommonPerformanceOptions(cli::AppPtr app,
                                 const ReadCollapserOptionsPtr& options,
                                 const std::string& group_name);

void AddConsensusOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddMarkdupOutputOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddConsensusOutputOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddConsensusClusterOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void AddConsensusDebugOptions(cli::AppPtr app, const ReadCollapserOptionsPtr& options, const std::string& group_name);

void ValidateConsensusOptions(const cli::ConstAppPtr& app, const ReadCollapserOptionsPtr& options);

void ValidateMarkdupOptions(const cli::ConstAppPtr& app, const ReadCollapserOptionsPtr& options);

}  // namespace xoos::read_collapser
