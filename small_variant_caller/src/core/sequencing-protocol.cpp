#include "sequencing-protocol.h"

namespace xoos::svc {

bool IsDuplexProtocol(const SequencingProtocol protocol) {
  return protocol == SequencingProtocol::kDuplex || protocol == SequencingProtocol::kDuplexSimplex;
}

}  // namespace xoos::svc
