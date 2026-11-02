#include "io/alignment.h"

#include <string>

#include <htslib/sam.h>

#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/read-name/read-name-parser.h>
#include <xoos/types/int.h>
#include <xoos/util/sequence-functions.h>
#include <xoos/yc-decode/yc-decoder.h>

namespace xoos::read_collapser {

Alignment::Alignment(io::Bam1Ptr record, const bool expect_umis, const bool ignore_read_name_parsing_errors)
    : record(std::move(record)),
      is_hairpin_read(bam_aux_get(this->record.get(), "YC") != nullptr),
      end_pos(static_cast<u32>(bam_endpos(this->record.get()))),
      original_flags(this->record->core.flag) {
  UpdateReadAdapterInfo(*this, expect_umis, ignore_read_name_parsing_errors);
}

Alignment::Alignment(io::Bam1Ptr record)
    : record(std::move(record)),
      is_hairpin_read(bam_aux_get(this->record.get(), "YC") != nullptr),
      end_pos(static_cast<u32>(bam_endpos(this->record.get()))),
      original_flags(this->record->core.flag) {
}

bool Alignment::IsPartial() const {
  return !is_5p_complete || !is_3p_complete;
}

bool Alignment::IsFivePrimeComplete() const {
  return is_5p_complete;
}

bool Alignment::IsThreePrimeComplete() const {
  return is_3p_complete;
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
    if ((original_flags & BAM_FREVERSE) != 0) {
      decoded_yc_tag.ReverseComplement(read_name);
    }
    base_types = decoded_yc_tag.GetBaseTypes();
  }
  return base_types;
}

void UpdateReadAdapterInfo(Alignment& alignment, const bool expect_umis, const bool ignore_parsing_errors) {
  const std::string_view qname{bam_get_qname(alignment.record.get())};
  // Parse the UMI from the read name. The UMI is expected to be at the end of the qname, separated by '|'.
  const read_name::AdapterInfo adapter_info = ignore_parsing_errors
                                                  ? read_name::ParseReadNameNoExcept(qname, expect_umis)
                                                  : read_name::ParseReadName(qname, expect_umis);
  // During alignment the direction of the read is determined, but the UMIs are not updated.
  // If the alignment is reverse, we swap the 5' and 3' UMIs. If the UMI are sequences,
  // we also need to reverse complement them.
  read_name::Umi umi5p = std::nullopt;
  read_name::Umi umi3p = std::nullopt;
  if (expect_umis) {
    umi5p = adapter_info.umi5p;
    umi3p = adapter_info.umi3p;
    if (alignment.IsReverse()) {
      umi5p = read_name::IsSequenceUmi(umi5p) ? sequence::ReverseComplement(umi5p.value()) : umi5p;
      umi3p = read_name::IsSequenceUmi(umi3p) ? sequence::ReverseComplement(umi3p.value()) : umi3p;
      std::swap(umi5p, umi3p);
    }
  }
  alignment.umi5p = umi5p;
  alignment.umi3p = umi3p;
  // If we expect UMIs, then the presence of UMI on the read ends determines whether the read is partial or full.
  // If we do not expect UMIs, then the presence of SID on the read ends determines whether the read is partial or full.
  alignment.is_5p_complete = expect_umis ? adapter_info.has_5p_umi : adapter_info.has_5p_sid;
  alignment.is_3p_complete = expect_umis ? adapter_info.has_3p_umi : adapter_info.has_3p_sid;
  if (alignment.IsReverse()) {
    std::swap(alignment.is_3p_complete, alignment.is_5p_complete);
  }
  // Special handling for hairpin reads:
  // All hairpin reads must have both SIDs because reads with truncated SIDs (part of the hairpin adapter)
  // are filtered out during demux. If not clustering by UMIs (SBX-D/FAST), it is safe to assume that all
  // hairpin reads are full.
  // For hairpin libraries with UMIs (I2X), we need to distinguish between reads with both UMIs and
  // reads with only one UMI in order to perform UMI clustering correctly. In this case, we do not
  // assume that all reads are full. Instead, we rely on the presence of UMIs to determine whether
  // the read is partial or full.
  if (alignment.is_hairpin_read && !expect_umis) {
    alignment.is_3p_complete = true;
    alignment.is_5p_complete = true;
  }
}

bool Alignment::HasSaTag() const {
  return bam_aux_get(record.get(), "SA") != nullptr;
}

}  // namespace xoos::read_collapser
