#pragma once

#include <compare>
#include <limits>
#include <string>
#include <string_view>

#include <fmt/format.h>

#include <xoos/types/int.h>

namespace xoos::genome {

/**
 * A genomic region using 0-based, half-open coordinates [start, end).
 *
 * This is the canonical internal coordinate convention for XOOS code:
 * BED files natively use it, and htslib's iterator coordinates (sam_itr_queryi)
 * use it as well. Samtools-style region *strings* ("chr:start-end") are instead
 * 1-based closed; use ParseSamtoolsLikeRegion to convert from that form.
 */
struct Region {
  // Chromosome / contig name.
  std::string chromosome;
  // 0-based, inclusive start position.
  s64 start;
  // 0-based, exclusive end position.
  s64 end;

  // Defaulted as a friend (non-member) so Region stays a plain aggregate while
  // the compiler generates the lexicographic comparisons (by chromosome, then
  // start, then end). A friend keeps SonarQube's "struct with no member methods"
  // rule satisfied without hand-writing branch-laden operators. Defaulting <=>
  // also implicitly declares a matching friend operator==, so equality is
  // available without spelling it out.
  friend std::strong_ordering operator<=>(const Region& lhs, const Region& rhs) = default;
};

// Length of the interval in bases (end - start).
inline s64 Length(const Region& region) {
  return region.end - region.start;
}

/**
 * Which samtools-style region forms a parse call is allowed to accept.
 *
 * All forms use 1-based closed coordinates on input; the parser converts them
 * to the 0-based half-open coordinates stored in Region. The mode only controls
 * which *shapes* of input are permitted, not how coordinates are interpreted.
 */
enum class RegionForm : u8 {
  // Only an explicit range "chr:start-end" is accepted. Whole-chromosome and
  // single-position inputs are rejected. This matches the historical behaviour
  // of the callers that required both endpoints.
  kRange,
  // "chr:start-end" or a single position "chr:start". Whole-chromosome inputs
  // (a bare name, or a name whose suffix is not coordinate-like) are rejected.
  kRangeOrPosition,
  // Any samtools form: whole chromosome "chr", single position "chr:start", or
  // range "chr:start-end".
  kAny,
};

/**
 * Parse a samtools-like region string into a 0-based half-open Region.
 *
 * Accepted forms (mirroring `samtools` region syntax), subject to `form`:
 *   - "chr"             whole chromosome  -> [0, INT64_MAX)   (kAny only)
 *   - "chr:start"       single position   -> [start-1, start) (kAny, kRangeOrPosition)
 *   - "chr:start-end"   inclusive range   -> [start-1, end)   (all forms)
 *
 * The numeric coordinates in the string are 1-based and the range is closed
 * (both endpoints included), exactly like samtools. They are converted to
 * 0-based half-open coordinates in the returned Region.
 *
 * Differences from samtools (intentional, stricter validation):
 *   - A start of 0 is rejected with a clear error. Samtools/htslib silently
 *     clamp a 1-based start of 0 to position 1; here it is treated as a bug in
 *     the caller's input, because a valid 1-based coordinate is never 0.
 *   - end < start is rejected.
 *   - Forms not permitted by `form` are rejected with a clear error, so callers
 *     that only ever expect an explicit range fail fast on malformed input
 *     instead of silently producing a whole-chromosome span.
 *
 * Commas in coordinates (e.g. "chr1:1,000-2,000") are tolerated, matching
 * samtools' lenient parsing.
 *
 * The single-argument overload accepts any samtools form (equivalent to passing
 * RegionForm::kAny); the two-argument overload restricts which forms are allowed.
 *
 * @param region region string in samtools format
 * @param form   which region forms to accept
 * @return Region with 0-based half-open coordinates
 * @throws std::runtime_error on malformed input, a start of 0, end < start, or
 *         a form not permitted by `form`
 */
Region ParseSamtoolsLikeRegion(const std::string& region);
Region ParseSamtoolsLikeRegion(const std::string& region, RegionForm form);

/**
 * Extract the contig (sequence) name from a samtools-like region string.
 *
 * Uses the same boundary resolution as ParseSamtoolsLikeRegion: the contig is
 * everything before the *last* colon, but only when the text after that colon
 * looks like a coordinate range. Otherwise the whole string is the contig name,
 * so sequence names containing colons (e.g. "HLA:A" or "HLA:A:1-100") are
 * handled consistently with the parser.
 *
 * Examples:
 *   "chr1"            -> "chr1"
 *   "chr1:100-200"    -> "chr1"
 *   "HLA:A"           -> "HLA:A"
 *   "HLA:A:1-100"     -> "HLA:A"
 *
 * @param region region string in samtools format
 * @return the contig name portion
 */
std::string ContigName(std::string_view region);

/**
 * Check whether the contig of a samtools-like region string equals a name.
 *
 * Equivalent to `ContigName(region) == contig` but without allocating, so it is
 * cheap enough to call once per region when scanning a whole file. The contig
 * boundary is resolved exactly as in ContigName and ParseSamtoolsLikeRegion, so
 * a name containing colons is matched as a whole.
 *
 * The match is exact, not a prefix match, and the coordinate suffix is optional.
 *
 * Examples:
 *   ("chrX:1-100", "chrX")   -> true
 *   ("chrX", "chrX")         -> true
 *   ("chrXY:1-100", "chrX")  -> false
 *   ("HLA:A:1-100", "HLA:A") -> true
 *
 * @param region region string in samtools format
 * @param contig contig name to compare against
 * @return true if the region's contig is exactly `contig`
 */
bool IsEqualContig(std::string_view region, std::string_view contig);

}  // namespace xoos::genome

/**
 * Format a Region back into a samtools-like region string.
 *
 * Renders the inverse of ParseSamtoolsLikeRegion: the stored 0-based half-open
 * [start, end) coordinates are emitted as 1-based closed "chr:start-end". A
 * whole-chromosome span (end == INT64_MAX, as produced for a bare "chr") is
 * rendered as just the contig name, so a round trip through the parser and
 * formatter preserves the original string for all three samtools forms.
 */
// Inherit the string_view formatter so its parse() handles the format spec; we
// only override how a Region is rendered. fmt mandates the method be named
// `format`, which violates the project's PascalCase rule (S100) for an external
// API contract we cannot rename.
template <>
struct fmt::formatter<xoos::genome::Region> : public fmt::formatter<std::string_view> {
  template <typename Context>
  auto format(const xoos::genome::Region& region, Context& ctx) const {  // NOSONAR(S100) fmt API contract
    if (region.end == std::numeric_limits<xoos::s64>::max()) {
      return fmt::formatter<std::string_view>::format(region.chromosome, ctx);
    }
    // 0-based inclusive start -> 1-based; 0-based exclusive end == 1-based inclusive end.
    return fmt::formatter<std::string_view>::format(
        fmt::format("{}:{}-{}", region.chromosome, region.start + 1, region.end), ctx);
  }
};
