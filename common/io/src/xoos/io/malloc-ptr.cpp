#include "xoos/io/malloc-ptr.h"

#include <cstdlib>
#include <gsl/pointers>

namespace xoos::io {

// S5008/M23_259: htslib allocates output buffers via malloc() and exposes them as void*;
// MallocDeleter must accept void* and use manual memory management at this C API boundary.
void MallocDeleter::operator()(const gsl::owner<void*> ptr) const {  // NOSONAR(S5008,M23_259)
  // This deleter is the RAII wrapper for malloc-allocated htslib buffers; calling
  // free here is its entire purpose.
  // NOLINTNEXTLINE(cppcoreguidelines-no-malloc,hicpp-no-malloc)
  free(ptr);  // NOSONAR(S1231,S984) — MallocDeleter is the RAII release for htslib malloc buffers
}

}  // namespace xoos::io
