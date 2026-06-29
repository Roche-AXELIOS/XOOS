#pragma once

#include <xoos/types/int.h>

#include "io/read-record.h"

namespace xoos::demux {

namespace metrics_constraints {
inline u32 max_sid_id_index = 0;
inline u32 max_logged_read_length = kMaxAlignLength;
};  // namespace metrics_constraints

}  // namespace xoos::demux
