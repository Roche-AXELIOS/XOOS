#pragma once

#include <optional>
#include <string>

#include <htslib/sam.h>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

namespace xoos::io {

// Fields extracted from a single @RG (read group) line in a BAM/SAM header.
// Any field absent from the header line is left as an empty string.
struct ReadGroup {
  // Read group identifier (ID tag). Required by the SAM spec for every @RG line.
  std::string id;
  // Sample name (SM tag).
  std::string sample_name;
  // Sequencing platform (PL tag).
  std::string platform;
  // Library identifier (LB tag).
  std::string library;
  // Platform model (PM tag). SBX instrument BAMs from ppa2bam carry the
  // sequencing run type here (e.g. "SBX-D") rather than in LB, as do some legacy
  // germline_wgs BAMs.
  std::string platform_model;
  // Platform unit (PU tag).
  std::string platform_unit;
};

// Parse every @RG line in the header into a ReadGroup, one entry per line, in
// header order. Missing tags (SM, PL, LB, PM, PU) are returned as empty strings.
// Returns an empty vector if the header has no @RG lines.
vec<ReadGroup> SamHdrParseReadGroups(sam_hdr_t* hdr);

// Parse the first @RG read group from a BAM file header.
// Returns std::nullopt if the header contains no @RG lines.
std::optional<ReadGroup> GetFirstReadGroupFromBam(const fs::path& bam_path);

}  // namespace xoos::io
