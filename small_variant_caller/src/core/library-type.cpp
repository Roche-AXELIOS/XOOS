#include "library-type.h"

#include <optional>

#include <xoos/enum/enum-util.h>
#include <xoos/error/error.h>
#include <xoos/io/library-type.h>

namespace xoos::svc {

namespace {
/**
 * @brief Translate a library type detected by common I/O code into the SVC domain enum.
 * @param library_type Library type detected from BAM metadata or read content.
 * @return Equivalent SVC library type.
 * @throws xoos::error::Error if the common I/O enum contains an unsupported value.
 */
LibraryType FromIoLibraryType(const io::LibraryType library_type) {
  using enum io::LibraryType;
  switch (library_type) {
    case kSbxS:
      return LibraryType::kSbxS;
    case kSbxD:
      return LibraryType::kSbxD;
    case kSbxFast:
      return LibraryType::kSbxFast;
    default:
      throw error::Error("Unknown BAM library type");
  }
}
}  // namespace

LibraryType ResolveLibraryType(const LibraryType requested, const vec<fs::path>& bam_paths) {
  // An explicit CLI value is authoritative. Avoid opening BAMs in this path so callers that do not have BAM input can
  // still select a profile directly.
  if (requested != LibraryType::kAutoDetect) {
    return requested;
  }
  if (bam_paths.empty()) {
    throw error::Error("Library type cannot be auto-detected without a BAM input; provide --library-type explicitly");
  }

  std::optional<LibraryType> resolved;
  // Multi-BAM commands must use one profile and model set, so every input must resolve to the same library type.
  for (const auto& bam_path : bam_paths) {
    const auto detected = io::DetectLibraryType(bam_path);
    if (!detected.has_value()) {
      throw error::Error("Could not automatically detect library type from BAM: {}", bam_path);
    }
    const auto current = FromIoLibraryType(*detected);
    if (resolved.has_value() && *resolved != current) {
      throw error::Error("BAM inputs have conflicting library types: '{}' and '{}'",
                         enum_util::FormatEnumName(*resolved),
                         enum_util::FormatEnumName(current));
    }
    resolved = current;
  }
  return *resolved;
}

std::string LibraryTypeModelToken(const LibraryType library_type) {
  using enum LibraryType;
  switch (library_type) {
    case kSbxS:
      return "sbxs";
    case kSbxD:
      return "sbxd";
    case kSbxFast:
      return "sbxfast";
    case kAutoDetect:
      throw error::Error("Cannot select a model while library type is auto-detect");
    default:
      throw error::Error("Unknown library type");
  }
}

}  // namespace xoos::svc
