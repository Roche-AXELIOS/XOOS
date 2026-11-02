#include "core/supplementary-alignment.h"

#include <cstddef>
#include <string>
#include <string_view>
#include <utility>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/types/vec.h>
#include <xoos/util/parse-int.h>
#include <xoos/util/string-functions.h>

namespace xoos::alignment_metrics {

namespace {

// The SA tag encodes each supplementary alignment as a comma-separated record:
//   rname,pos,strand,CIGAR,mapQ,NM
// These indices identify each field within a single record once split on ','.
constexpr size_t kFieldRname = 0;
constexpr size_t kFieldPos = 1;
constexpr size_t kFieldStrand = 2;
constexpr size_t kFieldCigar = 3;
constexpr size_t kFieldMapq = 4;
// Field 5 is the edit distance (NM). It is not stored, but is still required for
// the record to be well-formed (see kNumFields).
constexpr size_t kNumFields = 6;
// Placeholder rname for supplementary alignments that are unmapped or have missing reference names.
constexpr auto kUnmappedRname = "*";

}  // namespace

vec<SaTag> ParseSaTag(const std::string_view sa_tag_value) {
  vec<SaTag> sa_tags;
  // Each supplementary alignment is one ';'-delimited record. A trailing ';'
  // (the conventional SA tag terminator) yields an empty final entry that we skip.
  const auto records = string::Split(std::string(sa_tag_value), ";");
  for (const auto& record : records) {
    if (record.empty()) {
      continue;
    }
    const auto fields = string::Split(record, ",");
    if (fields.size() != kNumFields) {
      throw error::Error("Malformed SA tag entry '{}': expected {} comma-separated fields but found {}.",
                         record,
                         kNumFields,
                         fields.size());
    }

    const std::string& strand = fields[kFieldStrand];
    if (strand != "+" && strand != "-") {
      throw error::Error("Malformed SA tag entry '{}': strand field must be '+' or '-' but was '{}'.", record, strand);
    }

    SaTag sa_tag;
    sa_tag.rname = fields[kFieldRname];
    if (sa_tag.rname == kUnmappedRname) {
      continue;
    }
    // SA tag positions are 1-based; validate before converting to 0-based so a malformed pos
    // (0 or negative) cannot silently produce a negative coordinate.
    const s64 one_based_pos = util::ParseS64(fields[kFieldPos]);
    if (one_based_pos < 1) {
      throw error::Error("Malformed SA tag entry '{}': position field must be a 1-based coordinate >= 1 but was '{}'.",
                         record,
                         fields[kFieldPos]);
    }
    sa_tag.pos = one_based_pos - 1;
    sa_tag.is_reverse = strand == "-";
    sa_tag.cigar = io::SamParseCigar(fields[kFieldCigar]);
    sa_tag.mapq = util::ParseU8(fields[kFieldMapq]);
    sa_tags.push_back(std::move(sa_tag));
  }
  return sa_tags;
}

}  // namespace xoos::alignment_metrics
