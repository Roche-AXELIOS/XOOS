#include "xoos/read-name/umi.h"

#include <cstddef>
#include <optional>
#include <string>

#include "xoos/util/hash.h"

namespace xoos::read_name {

size_t UmiPairHash::operator()(const UmiPair& umi_pair) const {
  return util::hash::Hash(umi_pair.umi3p, umi_pair.umi5p);
}

bool IsSequenceUmi(const std::optional<std::string>& umi) {
  return umi.has_value() && umi->find_first_not_of("ACGT") == std::string::npos;
}

}  // namespace xoos::read_name
