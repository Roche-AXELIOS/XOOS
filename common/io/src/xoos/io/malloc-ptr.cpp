#include "xoos/io/malloc-ptr.h"

#include <cstdlib>

namespace xoos::io {

void MallocDeleter::operator()(gsl::owner<void*> ptr) {
  free(ptr);
}

}  // namespace xoos::io
