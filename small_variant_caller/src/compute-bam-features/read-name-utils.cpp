#include "read-name-utils.h"

#include <xoos/util/string-functions.h>

namespace xoos::svc {
void ParseReadName(const std::string& read_name, u32& plus_counts, u32& minus_counts, u32& family_size) {
  // Get cluster size and strand counts information
  // example read name: 1.1.1371c-3-4
  // 1.1.1371c is the read name
  // 3 is the number of reads supporting the variant on the minus strand
  // 4 is the total number of reads supporting the variant
  auto read_name_parts = string::Split(read_name, "-");
  family_size = std::stoi(*(read_name_parts.end() - 1));
  minus_counts = std::stoi(*(read_name_parts.end() - 2));
  plus_counts = family_size - minus_counts;
}

}  // namespace xoos::svc
