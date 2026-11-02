#include "xoos/io/bam/read-group.h"

#include <optional>
#include <string>

#include <htslib/sam.h>

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "xoos/io/htslib-util/htslib-util.h"
#include "xoos/io/htslib-util/kstring.h"

namespace xoos::io {

namespace {

// SAM @RG line type and its tags (see the SAM spec, section 1.3). Declared as
// null-terminated C strings so they pass straight to the htslib C API without
// an extra copy, while still avoiding array-to-pointer decay at call sites
// (the constants are already pointers, not arrays).
const constexpr char* kReadGroupLineType = "RG";
/// Read group identifier.
const constexpr char* kReadGroupIdTag = "ID";
/// Sample name.
const constexpr char* kReadGroupSampleTag = "SM";
/// Sequencing platform.
const constexpr char* kReadGroupPlatformTag = "PL";
/// Library identifier.
const constexpr char* kReadGroupLibraryTag = "LB";
/// Platform model.
const constexpr char* kReadGroupPlatformModelTag = "PM";
/// Platform unit.
const constexpr char* kReadGroupPlatformUnitTag = "PU";

/// Read a single tag value for the @RG line identified by `rg_id`.
/// Returns an empty string if the tag is absent (sam_hdr_find_tag_id returns
/// non-zero).
// S945: sam_hdr_find_tag_id requires NUL-terminated C strings; std::string_view does not
// guarantee NUL termination, making .data() unsafe at this boundary. const char* is correct here.
std::string ReadGroupTag(sam_hdr_t* const hdr,  // NOSONAR(S945)
                         const char* const rg_id,
                         const char* const tag) {
  Kstring value;
  if (sam_hdr_find_tag_id(hdr, kReadGroupLineType, kReadGroupIdTag, rg_id, tag, value.Get()) == 0) {
    return value.Str();
  }
  return {};
}

}  // namespace

vec<ReadGroup> SamHdrParseReadGroups(sam_hdr_t* const hdr) {
  vec<ReadGroup> read_groups;
  const auto line_count = sam_hdr_count_lines(hdr, kReadGroupLineType);
  if (line_count < 0) {
    throw error::Error("Failed to count @RG lines in BAM header");
  }
  for (s32 pos = 0; pos < line_count; ++pos) {
    // ID is required by the SAM spec and doubles as the key for tag lookups.
    // htslib rejects @RG lines without an ID at header-parse time, so this null
    // check is defensive against a malformed in-memory header.
    const char* const rg_id = sam_hdr_line_name(hdr, kReadGroupLineType, pos);
    if (rg_id == nullptr) {
      throw error::Error("@RG line {} in BAM header is missing the required ID tag", pos);
    }
    read_groups.push_back(ReadGroup{.id = rg_id,
                                    .sample_name = ReadGroupTag(hdr, rg_id, kReadGroupSampleTag),
                                    .platform = ReadGroupTag(hdr, rg_id, kReadGroupPlatformTag),
                                    .library = ReadGroupTag(hdr, rg_id, kReadGroupLibraryTag),
                                    .platform_model = ReadGroupTag(hdr, rg_id, kReadGroupPlatformModelTag),
                                    .platform_unit = ReadGroupTag(hdr, rg_id, kReadGroupPlatformUnitTag)});
  }
  return read_groups;
}

std::optional<ReadGroup> GetFirstReadGroupFromBam(const fs::path& bam_path) {
  const auto bam_file = HtsOpen(bam_path, "rb");
  const auto hdr = SamHdrRead(bam_file.get());
  auto read_groups = SamHdrParseReadGroups(hdr.get());
  if (read_groups.empty()) {
    return std::nullopt;
  }
  return read_groups.front();
}

}  // namespace xoos::io
