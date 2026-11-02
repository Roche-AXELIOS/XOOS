#pragma once

#include <functional>
#include <map>
#include <string>
#include <utility>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::cnc {

/**
 * @brief Position lookup for excluding putative SNPs that overlap blocklisted regions.
 *
 * Loads a BED file of 0-based half-open regions and answers overlap queries by
 * (chromosome, 0-based position). Intended to filter out putative het-SNP sites
 * before they are kept for BAF extraction.
 */
class SnpBlocklist {
 public:
  SnpBlocklist() = default;

  /// Load intervals from a BED file. Regions are merged per chromosome so overlap
  /// queries are correct even when the input contains overlapping or unsorted entries.
  explicit SnpBlocklist(const fs::path& bed_file);

  /// @return true if the 0-based position on the given chromosome falls within any
  /// blocklisted region (start inclusive, end exclusive).
  bool Contains(const std::string& chromosome, s64 position) const;

  /// @return true if no regions were loaded.
  bool Empty() const;

  /// @return total number of blocklisted regions across all chromosomes.
  size_t Size() const;

 private:
  // Per-chromosome intervals sorted by start, as (start, end) half-open pairs.
  // std::less<> enables heterogeneous lookup, avoiding a temporary std::string per query.
  std::map<std::string, vec<std::pair<s64, s64>>, std::less<>> _intervals;
};

}  // namespace xoos::cnc
