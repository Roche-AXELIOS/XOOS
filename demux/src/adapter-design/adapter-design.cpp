#include "adapter-design/adapter-design.h"

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/vfs/vfs-iterator.h>
#include <xoos/vfs/vfs.h>

#include <filesystem>
#include <format>
#include <memory>
#include <unordered_map>

#include <nlohmann/json.hpp>

namespace vfs = xoos::vfs;

using json = nlohmann::json;  // NOLINT

namespace xoos::demux {
// clang-format off
NLOHMANN_JSON_SERIALIZE_ENUM(
  BarcodeType,
  {
    {BarcodeType::kLoop, "loop"},
    {BarcodeType::kRunway, "runway"},
    {BarcodeType::kUmi, "umi"},
    {BarcodeType::kSid, "sid"},
    {BarcodeType::kStem, "stem"},
    {BarcodeType::kAnchor, "anchor"}
  }
)
// clang-format on

// clang-format off
NLOHMANN_JSON_SERIALIZE_ENUM(
  LutTransform,
  {
    {LutTransform::kNone, "none"},
    {LutTransform::kReverse, "reverse"},
    {LutTransform::kReverseComplement, "reverse_complement"},
  }
)
// clang-format on

// clang-format off
NLOHMANN_JSON_SERIALIZE_ENUM(AdapterType,
  {
    {AdapterType::kDuplex, "duplex"},
    {AdapterType::kDuplexUMI, "duplex_umi"},
    {AdapterType::kDuplexStem, "duplex_stem"},
    {AdapterType::kYsuSl, "YSU-SL"},
    {AdapterType::kYsuTe, "YSU-TE"},
    {AdapterType::kYs, "YS"}
  }
)
// clang-format on

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(LutDefinition, sequences, max_edit_distance, transform)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SearchDefinition, max_wiggle_left, max_wiggle_right)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(BarcodeDefinition, type, lut, search)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AdapterDesign, name, type, adapter_5p, adapter_3p)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(AdapterDesignManifest, default_adapter_design_name, adapter_designs)

std::string Format(BarcodeType type) {
  std::unordered_map<BarcodeType, std::string> barcode_type_names = {
      {BarcodeType::kRunway, "Runway"}, {BarcodeType::kSid, "Sid"},       {BarcodeType::kUmi, "Umi"},
      {BarcodeType::kStem, "Stem"},     {BarcodeType::kAnchor, "Anchor"}, {BarcodeType::kLoop, "Loop"}};
  return barcode_type_names.at(type);
}

std::string FormatJson(const BarcodeDefinition& bd) {
  json j;
  to_json(j, bd);
  return j.dump();
}

AdapterDesignManifest LoadBundleManifest(const vfs::VirtualFilesystemPtr& vfs, const fs::path& path);

AdapterDesign FindAdapterDesign(const AdapterDesignManifest& manifest, const std::string& design_name);

AdapterDesign LoadAdapterDesign(const fs::path& adapter_design_bundle,
                                const std::optional<std::string>& adapter_design_name) {
  vfs::VirtualFilesystemPtr vfs = vfs::Open(adapter_design_bundle);
  AdapterDesignManifest manifest = LoadBundleManifest(vfs, fs::path{"manifest.json"});
  std::string design_name = adapter_design_name.value_or(manifest.default_adapter_design_name);
  Logging::Info("Loading adapter design '{}' from bundle '{}'", design_name, adapter_design_bundle.string());
  return FindAdapterDesign(manifest, design_name);
}

/**
 * Load JSON manifest from bundle.
 */
AdapterDesignManifest LoadBundleManifest(const vfs::VirtualFilesystemPtr& vfs, const fs::path& path) {
  vfs::VirtualFileHandlePtr fh = vfs->Open(path);
  if (fh == nullptr) {
    throw error::Error("Manifest '{}' does not exist in bundle '{}'", path.string(), vfs->GetName());
  }

  AdapterDesignManifest manifest;
  try {
    from_json(json::parse(std::make_shared<vfs::VirtualFileStream>(fh)), manifest);
  } catch (const json::exception& ex) {
    throw error::Error("Manifest '{}' in bundle '{}' cannot be parsed due to '{}'", path.string(), vfs->GetName(),
                       ex.what());
  }
  return manifest;
}

/**
 * Find adapter design with provided name, if there are multiple adapter design with the same name
 * only the first will be found.
 */
AdapterDesign FindAdapterDesign(const AdapterDesignManifest& manifest, const std::string& design_name) {
  auto adapter_design = std::ranges::find_if(manifest.adapter_designs,
                                             [&design_name](const auto& item) { return item.name == design_name; });
  if (adapter_design == manifest.adapter_designs.end()) {
    std::string available_designs;
    for (const auto& design : manifest.adapter_designs) {
      available_designs += std::format(" {} type: {}", design.name, static_cast<u32>(design.type));
    }
    throw error::Error("Adapter design '{}' does not exist in manifest, available designs:{}", design_name,
                       available_designs);
  }
  return *adapter_design;
}
}  // namespace xoos::demux
