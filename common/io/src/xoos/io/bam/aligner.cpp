#include "xoos/io/bam/aligner.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <locale>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/io/htslib-util/kstring.h>
#include <xoos/types/int.h>
#include <xoos/util/string-functions.h>

namespace xoos::io {

namespace {

// @PG tokens that identify a known aligner family. Parabricks records the invoked base aligner as a `pbrun`
// subcommand in the PN and ID tags (observed on real BAMs as `PN:pbrun fq2bam` and `PN:pbrun giraffe`), so the
// subcommand tokens (`fq2bam`, `giraffe`, `minimap2`) are matched in addition to the standalone aligner names.
// `vg` produces giraffe alignments (giraffe records `ID:0 PN:vg`, so the PN tag must be read); `bwa`/`bwa-mem2`/
// `fq2bam` produce bwa-style linear alignments; `minimap2` produces minimap2 alignments.
constexpr auto* kProgramGiraffe = "giraffe";
constexpr auto* kProgramVg = "vg";
constexpr auto* kProgramBwa = "bwa";
constexpr auto* kProgramFq2Bam = "fq2bam";
constexpr auto* kProgramMinimap2 = "minimap2";

std::string ToLower(std::string value) {
  std::ranges::transform(value, value.begin(), [](const char c) { return std::tolower(c, std::locale::classic()); });
  return value;
}

// Map a @PG field value to a known aligner family. Matching is case-insensitive and substring-based to tolerate
// variations in how aligners record their program name, and to match parabricks `pbrun <subcommand>` invocations
// recorded in the command line. Returns kUnknown when no aligner token is present.
Aligner AlignerFromProgramField(const std::string& value) {
  using enum Aligner;
  const std::string lowered = ToLower(value);
  if (lowered.contains(kProgramGiraffe) || lowered == kProgramVg) {
    return kGiraffe;
  }
  // minimap2 must be checked before bwa: parabricks fq2bam maps to bwa, but a pbrun minimap2 command line does
  // not contain a bwa token, so ordering here only guards against unrelated substrings.
  if (lowered.contains(kProgramMinimap2)) {
    return kMinimap2;
  }
  if (lowered.contains(kProgramBwa) || lowered.contains(kProgramFq2Bam)) {
    return kBwa;
  }
  return kUnknown;
}

// Read a tag value from a specific @PG record, returning an empty string when the tag is absent.
std::string ProgramTag(sam_hdr_t* const header, const char* const id, const char* const tag) {
  Kstring tag_value;
  if (sam_hdr_find_tag_id(header, io::kPgTag, io::kPgProgramIdTag, id, tag, tag_value.Get()) == 0) {
    return tag_value.Str();
  }
  return {};
}

// Number of leading space-delimited tokens of the command line inspected for an aligner token. The aligner
// invocation always appears at the start of the command line (e.g. `bwa mem ...`, `minimap2 -a ...`,
// `pbrun giraffe ...`), so only the program and its immediate subcommand are considered. The remainder of the
// command line (paths, option values, reference/index filenames) is ignored to avoid false-positive matches.
constexpr std::size_t kCommandLineTokenLimit = 2;

// Return the first `kCommandLineTokenLimit` space-delimited tokens of `value`, joined by single spaces.
std::string LeadingCommandTokens(const std::string& value) {
  const std::vector<std::string> tokens = string::Split(value, " ");
  const std::size_t count = std::min(tokens.size(), kCommandLineTokenLimit);
  return string::Join({tokens.begin(), tokens.begin() + static_cast<std::ptrdiff_t>(count)}, " ");
}

}  // namespace

Aligner DetectAlignerFromHeader(sam_hdr_t* const header) {
  using enum Aligner;
  s32 pos = 0;
  while (const char* const id = sam_hdr_line_name(header, io::kPgTag, pos)) {
    // Check the program name (PN) tag, then the program record ID, then the command line (CL). Parabricks
    // records the base aligner subcommand in PN/ID (`pbrun fq2bam`, `pbrun giraffe`), which the PN check
    // catches; CL is inspected as a defensive fallback for tools that only record the aligner there. Only the
    // leading tokens of the CL are inspected so that arguments (paths, ref/index filenames) cannot produce a
    // false-positive aligner match, whereas PN/ID are short identifiers matched in full.
    const std::array<std::string, 3> fields{ProgramTag(header, id, io::kPgProgramNameTag),
                                            id,
                                            LeadingCommandTokens(ProgramTag(header, id, io::kPgCommandLineTag))};
    for (const std::string& field : fields) {
      if (const Aligner aligner = AlignerFromProgramField(field); aligner != kUnknown) {
        return aligner;
      }
    }
    ++pos;
  }
  return kUnknown;
}

}  // namespace xoos::io
