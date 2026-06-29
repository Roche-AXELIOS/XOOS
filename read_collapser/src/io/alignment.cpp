#include "io/alignment.h"

#include <string>

#include <htslib/sam.h>

#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/types/int.h>
#include <xoos/yc-decode/yc-decoder.h>

namespace xoos::read_collapser {

Alignment::Alignment(io::Bam1Ptr record)
    : record(std::move(record)),
      end_pos(static_cast<u32>(bam_endpos(this->record.get()))),
      original_flags(this->record->core.flag) {
}

bool Alignment::IsPartial() const {
  return (umi5p == std::nullopt || umi3p == std::nullopt) && !(umi5p == std::nullopt && umi3p == std::nullopt);
}

bool Alignment::IsReverse() const {
  return bam_is_rev(record.get());
}

bool Alignment::IsForward() const {
  return !IsReverse();
}

u32 Alignment::StartPos() const {
  return static_cast<u32>(record->core.pos);
}

u32 Alignment::EndPos() const {
  return end_pos;
}

vec<yc_decode::BaseType> Alignment::GetBaseTypes() const {
  auto yc_tag = io::BamAuxGet<std::string>(this->record.get(), "YC");
  vec<yc_decode::BaseType> base_types;
  if (yc_tag != std::nullopt) {
    const std::string read_name = bam_get_qname(this->record.get());
    // We pass the YC tag string instead of the BAM record because yc_decode would try to reverse
    // the YC tag if it sees the BAM record has the reverse strand flag set. We want to avoid this
    // behavior because we manually set the strand for parent-parent duplex during deconvolution and
    // this behavior would mess up the order of the base types for R2 reads in the parent-parent workflow.
    auto decoded_yc_tag = yc_decode::DeserializeYcTag(*yc_tag, read_name);
    // If the original BAM record is on the reverse strand, we need to reverse complement the decoded YC tag
    // to get the correct order of base types
    if (original_flags & BAM_FREVERSE) {
      decoded_yc_tag.ReverseComplement(read_name);
    }
    base_types = decoded_yc_tag.GetBaseTypes();
  }
  return base_types;
}

}  // namespace xoos::read_collapser
