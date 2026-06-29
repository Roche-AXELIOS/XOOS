#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace xoos::demux {
enum class BarcodeType {
  kRunway,
  kUmi,
  kSid,
  kStem,
  kAnchor,
  kLoop,
};

std::string Format(BarcodeType type);

enum class LutTransform {
  kReverse,
  kReverseComplement,
  kNone,
};

/**
 * Define how the LUT is generated from a set of input sequences.
 */
struct LutDefinition {
  fs::path sequences;        /// A FASTA file containing all of the sequences and their names for this LUT
  int max_edit_distance{0};  /// The LUT will contain all sequences within this max_edit_distance of the sequences
  LutTransform transform{LutTransform::kNone};  /// The sequences can first be transformed before generating the LUT

  friend bool operator==(const LutDefinition& a, const LutDefinition& b) = default;
};

/**
 * Define how to search for a barcode in an adapter, search to the left and right of the expected position by the
 * defined amount.
 */
struct SearchDefinition {
  int max_wiggle_left;
  int max_wiggle_right;

  friend bool operator==(const SearchDefinition& a, const SearchDefinition& b) = default;
};

struct BarcodeDefinition {
  BarcodeType type{};
  LutDefinition lut;
  SearchDefinition search{};

  friend bool operator==(const BarcodeDefinition& a, const BarcodeDefinition& b) = default;
};

std::string FormatJson(const BarcodeDefinition& bd);

enum class AdapterType { kDuplex, kDuplexUMI, kDuplexStem, kYsuSl, kYsuTe, kYs };

struct AdapterDesign {
  std::string name;
  AdapterType type{};
  std::vector<BarcodeDefinition> adapter_5p;
  std::vector<BarcodeDefinition> adapter_3p;
};

/**
 * A struct used to describe multiple potential adapter designs.
 * This struct will most often be loaded from a bundle containing multiple adapter designs,
 * it will describe the details of the demultiplexing and trimming logic used for that adapter design.
 */
struct AdapterDesignManifest {
  std::string default_adapter_design_name;
  std::vector<AdapterDesign> adapter_designs;
};

/**
 * Load adapter design with the given name from the given bundle, or use default if name is nullopt.
 *
 * Bundle must contain a "manifest.json" which contains an AdapterDesignManifest in JSON representation.
 */
AdapterDesign LoadAdapterDesign(const fs::path& adapter_design_bundle,
                                const std::optional<std::string>& adapter_design_name);
}  // namespace xoos::demux
