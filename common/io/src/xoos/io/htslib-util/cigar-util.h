#pragma once

#include <xoos/types/int.h>

namespace xoos::io {

// Whether the CIGAR op consumes query/read bases (M, I, S, =, X).
bool ConsumesRead(u32 op);  // NOSONAR(S5536)

// Whether the CIGAR op consumes reference bases (M, D, N, =, X).
bool ConsumesReference(u32 op);  // NOSONAR(S5536)

}  // namespace xoos::io
