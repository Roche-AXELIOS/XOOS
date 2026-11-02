#include "xoos/io/library-type.h"

#include <algorithm>
#include <array>
#include <gsl/gsl>
#include <locale>
#include <optional>
#include <set>
#include <string>
#include <string_view>

#include <fmt/format.h>
#include <htslib/sam.h>

#include <xoos/enum/enum-util.h>
#include <xoos/io/bam/read-group.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/float.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/util/string-functions.h>
#include <xoos/yc-decode/yc-decoder.h>

namespace xoos::io {

namespace {

/// Raw Phred values represented by '+', '/', '0', and '5' in Phred+33. These
/// values come from the current Everest Output Spec read specification and are
/// subject to change.
constexpr std::array<u8, 4> kMultiValueSimplexBaseQualities{10, 14, 15, 20};

// The LB and PU formats follow the SBX SAM Format Spec's Read Group Header
// Design. The separator and adapter-class indexes may change in future revisions.
constexpr char kReadGroupFieldSeparator = '.';
constexpr size_t kPlatformUnitAdapterClassIndex = 4;
constexpr size_t kLibraryAdapterClassIndex = 1;

/// Result of scanning the first ::kLibraryTypeReadSampleSize reads.
struct BamReadPreviewResults {
  /// At least one read was sampled.
  bool saw_reads{false};
  /// At least one sampled read carried a YC tag.
  bool saw_yc_tag{false};
  /// At least one YC tag encoded a duplex read.
  bool saw_duplex_yc{false};
  /// Total bases seen across sampled reads.
  u64 base_count{0};
  /// Bases with BQ > kLibraryTypeHighBqThreshold.
  u64 high_bq_base_count{0};
  /// Distinct base qualities across sampled reads with quality strings.
  std::set<u8> observed_base_qualities;

  /// Fraction of sampled bases exceeding the high-BQ threshold.
  f64 HighBqFraction() const {
    return base_count == 0 ? 0.0 : static_cast<f64>(high_bq_base_count) / static_cast<f64>(base_count);
  }

  /// Whether the sampled qualities match a supported simplex profile.
  bool HasSimplexBaseQualities() const {
    // XOOS Demux can override simplex base qualities with a single configured
    // value, and older SBX-S data used uniform qualities such as Q22 or Q18.
    if (observed_base_qualities.size() == 1) {
      return true;
    }
    return !observed_base_qualities.empty() &&
           observed_base_qualities.size() <= kMultiValueSimplexBaseQualities.size() &&
           std::ranges::all_of(observed_base_qualities, [](const u8 qual) {
             return std::ranges::find(kMultiValueSimplexBaseQualities, qual) != kMultiValueSimplexBaseQualities.end();
           });
  }
};

/// Upper-case a copy of the string for case-insensitive run_type matching.
std::string ToUpper(std::string value) {
  std::ranges::transform(
      value, value.begin(), [](const char chr) { return std::toupper(chr, std::locale::classic()); });
  return value;
}

/// Human-readable library type name for logging (e.g. "sbx-d").
std::string LibraryTypeName(const LibraryType library_type) {
  return enum_util::FormatEnumName(library_type);
}

/// Extract and convert the adapter-class component of a structured @RG field.
std::optional<LibraryType> LibraryTypeFromStructuredField(const std::string_view value,
                                                          const size_t adapter_class_index) {
  const auto components = string::SplitStringView(value, kReadGroupFieldSeparator);
  if (components.size() <= adapter_class_index) {
    return std::nullopt;
  }
  return LibraryTypeFromRunType(std::string{components[adapter_class_index]});
}

/// Scan the first ::kLibraryTypeReadSampleSize reads of a BAM, collecting the
/// content signals used by the detection fallbacks and the consistency check.
BamReadPreviewResults ScanReadContent(const fs::path& bam_path) {
  BamReadPreviewResults signals;
  // "r" auto-detects the format so SAM/BAM/CRAM inputs all read correctly.
  const auto bam_file = SamOpen(bam_path, "r");
  const auto hdr = SamHdrRead(bam_file.get());
  const auto record = BamInit1();

  u32 read_count = 0;
  while (read_count < kLibraryTypeReadSampleSize) {
    if (!SamRead1(bam_file.get(), hdr.get(), record.get())) {
      break;
    }
    ++read_count;
    signals.saw_reads = true;

    // YC tag: presence and duplex encoding. DeserializeYcTag returns an empty
    // (simplex) tag when the record carries no YC aux field.
    if (bam_aux_get(record.get(), "YC") != nullptr) {
      signals.saw_yc_tag = true;
      const auto yc_tag = yc_decode::DeserializeYcTag(record.get());
      if (yc_tag.IsDuplex()) {
        signals.saw_duplex_yc = true;
      }
    }

    // Base-quality distribution. htslib stores 0xff in the first byte when the
    // record has no quality string ("*"); such records are skipped.
    static constexpr u8 kMissingQual = 0xff;
    const s32 len = record->core.l_qseq;
    const u8* const quals_ptr = bam_get_qual(record.get());
    if (len > 0 && quals_ptr != nullptr && *quals_ptr != kMissingQual) {
      const gsl::span<const u8> quals{quals_ptr, static_cast<size_t>(len)};
      for (const u8 qual : quals) {
        signals.observed_base_qualities.insert(qual);
        ++signals.base_count;
        if (qual > kLibraryTypeHighBqThreshold) {
          ++signals.high_bq_base_count;
        }
      }
    }
  }
  return signals;
}

}  // namespace

std::optional<LibraryType> LibraryTypeFromRunType(const std::string& run_type) {
  using enum LibraryType;
  // The pipeline stamps run_type verbatim into @RG LB (see xoosnf
  // subworkflows/local/align/nextflow.config). Known run_type values collapse
  // onto three behavioural representatives; versioned I2X designs match their
  // stable family prefix. Values are upper-cased for case-insensitive matching.
  const std::string normalized_run_type = ToUpper(run_type);
  if (normalized_run_type.starts_with(kRunTypeI2xHpsidPrefix) ||
      normalized_run_type.starts_with(kRunTypeI2xHpuPrefix)) {
    return kSbxD;
  }

  // Keyed on the canonical non-I2X run_type constants, upper-cased at
  // construction so lookups can match case-insensitively via ToUpper.
  static const StrUnorderedMap<LibraryType> kRunTypeToLibraryType = {
      // SBX-FAST is its own representative.
      {ToUpper(std::string{kRunTypeSbxFast}), kSbxFast},
      // Other duplex designs collapse to kSbxD.
      {ToUpper(std::string{kRunTypeSbxD}), kSbxD},
      {ToUpper(std::string{kRunTypeSbxDm}), kSbxD},
      {ToUpper(std::string{kRunTypeI2xI7HpsidDuplexV0}), kSbxD},
      // Simplex designs collapse to kSbxS.
      {ToUpper(std::string{kRunTypeYs}), kSbxS},
      {ToUpper(std::string{kRunTypeYsu}), kSbxS},
      {ToUpper(std::string{kRunTypeYsNew}), kSbxS},
  };
  const auto match = kRunTypeToLibraryType.find(normalized_run_type);
  if (match == kRunTypeToLibraryType.end()) {
    return std::nullopt;
  }
  return match->second;
}

bool IsDuplex(const LibraryType library_type) {
  using enum LibraryType;
  switch (library_type) {
    case kSbxD:
    case kSbxFast:
      return true;
    case kSbxS:
    // The simplex representative and the default share this branch: the default
    // guards an out-of-range value cast into the enum without a cloned return
    // that bugprone-branch-clone would flag.
    default:
      return false;
  }
}

std::optional<LibraryType> DetectLibraryType(const fs::path& bam_path) {
  // Steps 1 through 4: header read groups. The current output specification
  // embeds the adapter class in PU and LB; older inputs may store the run type
  // directly in LB or PM. Field priority is structured LB, full LB, structured
  // PU, then full PM.
  {
    // "r" auto-detects the format so SAM/BAM/CRAM inputs all read correctly.
    const auto bam_file = SamOpen(bam_path, "r");
    const auto hdr = SamHdrRead(bam_file.get());
    const auto read_groups = SamHdrParseReadGroups(hdr.get());

    // Step 1: structured LB adapter class.
    for (const auto& read_group : read_groups) {
      if (const auto library_type = LibraryTypeFromStructuredField(read_group.library, kLibraryAdapterClassIndex)) {
        Logging::Info("Detected library type {} from adapter-class component @RG LB[{}] of '{}'",
                      LibraryTypeName(*library_type),
                      kLibraryAdapterClassIndex,
                      bam_path.string());
        return library_type;
      }
    }

    // Step 2: full LB value.
    for (const auto& read_group : read_groups) {
      if (const auto library_type = LibraryTypeFromRunType(read_group.library)) {
        Logging::Info(
            "Detected library type {} from full @RG LB tag of '{}'", LibraryTypeName(*library_type), bam_path.string());
        return library_type;
      }
    }

    // Step 3: structured PU adapter class.
    for (const auto& read_group : read_groups) {
      if (const auto library_type =
              LibraryTypeFromStructuredField(read_group.platform_unit, kPlatformUnitAdapterClassIndex)) {
        Logging::Info("Detected library type {} from adapter-class component @RG PU[{}] of '{}'",
                      LibraryTypeName(*library_type),
                      kPlatformUnitAdapterClassIndex,
                      bam_path.string());
        return library_type;
      }
    }

    // Step 4: full PM value from instrument or legacy headers.
    for (const auto& read_group : read_groups) {
      if (const auto library_type = LibraryTypeFromRunType(read_group.platform_model)) {
        Logging::Warn(
            "Detected library type {} from full @RG PM tag (instrument/legacy location) of '{}'; PU and LB did not "
            "resolve the library type",
            LibraryTypeName(*library_type),
            bam_path.string());
        return library_type;
      }
    }
  }

  // Steps 5 through 7: read-content fallbacks over the first N reads.
  const BamReadPreviewResults signals = ScanReadContent(bam_path);

  // Steps 5 and 7 can only establish that the data is duplex, not which duplex
  // library type it is: SBX-D and SBX-FAST are indistinguishable from read
  // content alone. We return SBX-D as the representative duplex library type and
  // say so in the log so callers do not read too much into the specific value.

  // Step 5: any duplex-encoded YC tag implies a duplex library type.
  if (signals.saw_duplex_yc) {
    Logging::Warn(
        "Duplex library type assigned due to detection of YC tags in the first {} reads of '{}'.  Assigning SBX-D, but "
        "may also be SBX-FAST; @RG LB/PM "
        "tags did not resolve the library type",
        kLibraryTypeReadSampleSize,
        bam_path.string());
    return LibraryType::kSbxD;
  }

  // Step 6: a supported simplex BQ profile across reads without YC tags implies a simplex library type.
  if (!signals.saw_yc_tag && signals.HasSimplexBaseQualities()) {
    Logging::Warn(
        "Detected simplex library type from the base-quality distribution of '{}': sampled bases contain {} distinct "
        "simplex-compatible BQ value(s) and no sampled reads carry YC tags; @RG tags did not resolve the library type. "
        "Returning SBX-S",
        bam_path.string(),
        signals.observed_base_qualities.size());
    return LibraryType::kSbxS;
  }

  // Step 7: a high fraction of high-BQ bases implies the data is not simplex.
  if (signals.base_count > 0 && signals.HighBqFraction() > kLibraryTypeHighBqFraction) {
    Logging::Warn(
        "Detected duplex library type from the base-quality distribution of '{}': {:.2f}% of sampled bases exceed BQ "
        "{} (> {:.0f}% threshold); @RG tags and YC tags did not resolve the library type. Returning SBX-D as the "
        "representative duplex library type (may also be SBX-FAST)",
        bam_path.string(),
        signals.HighBqFraction() * 100.0,
        kLibraryTypeHighBqThreshold,
        kLibraryTypeHighBqFraction * 100.0);
    return LibraryType::kSbxD;
  }

  Logging::Warn(
      "Unable to detect library type for '{}': no usable @RG PU/LB/PM value, no duplex YC tags, no YC-free supported "
      "BQ profile, and high-base-quality distribution was inconclusive",
      bam_path.string());
  return std::nullopt;
}

std::optional<bool> IsDuplex(const fs::path& bam_path) {
  const auto library_type = DetectLibraryType(bam_path);
  if (!library_type) {
    return std::nullopt;
  }
  return std::optional<bool>{IsDuplex(*library_type)};
}

std::optional<std::string> CheckLibraryTypeConsistency(const fs::path& bam_path, const LibraryType declared) {
  const BamReadPreviewResults signals = ScanReadContent(bam_path);
  if (!signals.saw_reads) {
    // No reads to cross-check against; nothing to contradict.
    return std::nullopt;
  }

  if (IsDuplex(declared)) {
    // Declared duplex but no YC tags in the sampled reads.
    if (!signals.saw_yc_tag) {
      auto message = fmt::format(
          "Library type mismatch in '{}': header declares duplex ({}) but the first {} reads carry no YC tags",
          bam_path.string(),
          LibraryTypeName(declared),
          kLibraryTypeReadSampleSize);
      Logging::Error("{}", message);
      return message;
    }
  } else {
    // Declared simplex but the reads show duplex signals.
    if (signals.saw_duplex_yc) {
      auto message = fmt::format(
          "Library type mismatch in '{}': header declares simplex ({}) but the first {} reads carry duplex-encoded YC "
          "tags",
          bam_path.string(),
          LibraryTypeName(declared),
          kLibraryTypeReadSampleSize);
      Logging::Error("{}", message);
      return message;
    }
    if (signals.base_count > 0 && signals.HighBqFraction() > kLibraryTypeHighBqFraction) {
      auto message = fmt::format(
          "Library type mismatch in '{}': header declares simplex ({}) but {:.2f}% of sampled bases exceed BQ {} "
          "(> {:.0f}% threshold), indicating duplex data",
          bam_path.string(),
          LibraryTypeName(declared),
          signals.HighBqFraction() * 100.0,
          kLibraryTypeHighBqThreshold,
          kLibraryTypeHighBqFraction * 100.0);
      Logging::Error("{}", message);
      return message;
    }
  }

  return std::nullopt;
}

}  // namespace xoos::io
