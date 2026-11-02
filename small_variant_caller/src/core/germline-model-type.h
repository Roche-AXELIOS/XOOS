#pragma once

#include <xoos/types/int.h>

namespace xoos::svc {

/**
 * @brief Bundled germline WGS model family used for automatic model selection.
 * @details This enum controls which pre-trained SNV/indel model pair is selected when the user does not provide
 * explicit model paths. It is intentionally separate from LibraryType: library type selects the assay-specific feature
 * profile and filename token, while model type selects between the large and small bundled model families.
 */
enum class GermlineModelType : u8 {
  kLarge,
  kSmall,
};

}  // namespace xoos::svc
