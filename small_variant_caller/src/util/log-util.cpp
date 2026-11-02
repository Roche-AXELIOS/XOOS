#include "log-util.h"

namespace xoos::svc {

/**
 * @brief Modifiable flag to control whether to treat warn messages as errors.
 * If set to true, any warn messages will throw a runtime error.
 */
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables) — intentionally mutable runtime flag, set from CLI
bool warn_as_error = false;  // NOSONAR(cpp:S5421,cpp:M23_388) — intentionally mutable runtime flag, set from CLI

}  // namespace xoos::svc
