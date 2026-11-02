#pragma once

#include <string>
#include <string_view>

#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::alignment_metrics {

/**
 * @brief A single supplementary alignment parsed from a read's SA tag.
 *
 * The SA tag encodes one or more supplementary alignments for a read as a
 * semicolon-delimited list, where each entry has the form:
 *   rname,pos,strand,CIGAR,mapQ,NM;
 *
 * The fields used by this module are stored after parsing (the NM edit distance is
 * not retained). The CIGAR is decoded into htslib-encoded operations (the same
 * packed `op|len << 4` representation returned by `bam_get_cigar`) so it can be
 * inspected with the standard `bam_cigar_*` macros.
 *
 * The query coordinates implied by the CIGAR are expressed in the orientation of
 * the supplementary alignment itself (i.e. relative to `is_reverse`), not in the
 * orientation of the primary alignment.
 */
struct SaTag {
  // Reference/contig name the supplementary alignment maps to.
  std::string rname;
  // 0-based leftmost reference position of the supplementary alignment.
  s64 pos{};
  // True if the supplementary alignment is on the reverse strand.
  bool is_reverse{};
  // htslib-encoded CIGAR operations for the supplementary alignment.
  vec<u32> cigar;
  // Mapping quality of the supplementary alignment.
  u8 mapq{};
};

/**
 * @brief Criteria for excluding supplementary alignments from soft-clip coverage.
 *
 * A supplementary alignment described by an SA tag is only counted when its mapping
 * quality is at least `min_mapq` and none of its BAM flags overlap `exclude_flags`.
 */
struct SaTagFilter {
  // Minimum mapping quality a supplementary alignment must meet to be counted.
  u32 min_mapq{0};
  // Bitmask of BAM flags that disqualify a supplementary alignment when any are set.
  u32 exclude_flags{0};
};

/**
 * @brief Parse the value of a read's SA tag into a list of supplementary alignments.
 *
 * @param sa_tag_value The raw SA tag string, e.g. "chr1,100,+,50M50S,60,2;chr2,200,-,40S60M,55,1;".
 *                      A trailing semicolon is optional and empty entries are ignored.
 * @return One SaTag per supplementary alignment described by the tag.
 * @throws xoos::error::Error if an entry does not contain the expected six fields
 *         or if a numeric/CIGAR field cannot be parsed.
 */
vec<SaTag> ParseSaTag(std::string_view sa_tag_value);

}  // namespace xoos::alignment_metrics
