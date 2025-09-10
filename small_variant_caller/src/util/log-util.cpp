#include "log-util.h"

namespace xoos::svc {

/**
 * @brief Modifiable flag to control whether to treat warn messages as errors.
 * If set to true, any warn messages will throw a runtime error.
 */
bool warn_as_error = false;

}  // namespace xoos::svc
