#pragma once

#include <string>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::svc {

/**
 * @brief Library preparation type used to select germline WGS configuration and models.
 * @details This SVC-owned enum keeps CLI state, including auto-detection, separate from the lower-level library type
 * returned by BAM inspection in common/io.
 */
enum class LibraryType : u8 {
  kAutoDetect,
  kSbxS,
  kSbxD,
  kSbxFast,
};

/**
 * @brief Resolve an explicitly selected or auto-detected library type for one or more BAMs.
 * @details A concrete request is authoritative and bypasses BAM inspection. Auto-detection inspects every BAM and
 * requires all detected library types to agree.
 * @param requested Library type requested by the user.
 * @param bam_paths BAM paths used only when @p requested is LibraryType::kAutoDetect.
 * @return The concrete requested or detected library type.
 * @throws xoos::error::Error if auto-detection has no BAM, cannot identify a BAM, or finds conflicting types.
 */
LibraryType ResolveLibraryType(LibraryType requested, const vec<fs::path>& bam_paths);

/**
 * @brief Render the legacy filename token used by bundled germline models.
 * @details Canonical CLI values contain hyphens, while existing model filenames use compact tokens such as `sbxd`.
 * @param library_type Concrete library type to render.
 * @return Filename token for @p library_type.
 * @throws xoos::error::Error if @p library_type is auto-detect or is not recognized.
 */
std::string LibraryTypeModelToken(LibraryType library_type);

}  // namespace xoos::svc
