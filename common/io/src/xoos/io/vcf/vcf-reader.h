#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <xoos/io/htslib-util/htslib-bcf-ptr.h>
#include <xoos/io/htslib-util/kstring.h>
#include <xoos/types/str-container.h>

#include "vcf-header.h"
#include "vcf-record.h"

namespace xoos::io {

namespace fs = std::filesystem;

class VcfReader {
 public:
  explicit VcfReader(const fs::path& input_vcf_name);

  VcfHeaderPtr GetHeader() const {
    return _hdr;
  }

  StrMap<s32> GetContigIndexes() const;
  VcfRecordPtr GetNextRecord(s32 max_unpack) const;
  VcfRecordPtr GetNextRecord() const;
  std::vector<VcfRecordPtr> GetAllRecords() const;
  /// Set the query region using 0-based half-open coordinates [start, end).
  bool SetRegion(const std::string& chrom, u64 start, u64 end);
  /// Set the query region by target ID using 0-based half-open coordinates [start, end).
  bool SetRegion(s32 tid, u64 start, u64 end);
  VcfRecordPtr GetNextRegionRecord(s32 unpack);
  bool HasIndex();

 private:
  VcfHeaderPtr _hdr;
  HtsFileSharedPtr _input_file;
  HtsIdxPtr _hts_idx;
  TbxPtr _tbx_idx;
  Kstring _kstr;
  HtsItrPtr _hts_itr;
};

}  // namespace xoos::io
