#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <xoos/io/htslib-util/kstring.h>

#include "vcf-header.h"
#include "vcf-record.h"

namespace xoos::io {

struct TbxDeleter {
  void operator()(tbx_t* tbx) const;
};

using TbxPtr = std::unique_ptr<tbx_t, TbxDeleter>;

class VcfReader {
 public:
  explicit VcfReader(const std::string& input_vcf_name);

  VcfHeaderPtr GetHeader() const {
    return _hdr;
  }

  std::map<std::string, s32> GetContigIndexes();
  VcfRecordPtr GetNextRecord(int max_unpack) const;
  VcfRecordPtr GetNextRecord() const;
  std::vector<VcfRecordPtr> GetAllRecords() const;
  bool SetRegion(const std::string& chrom, u64 start, u64 end);
  bool SetRegion(int tid, u64 start, u64 end);
  VcfRecordPtr GetNextRegionRecord(int unpack);
  bool HasIndex();

 private:
  VcfHeaderPtr _hdr;
  HtsFileSharedPtr _input_file;
  HtsIdxPtr _hts_idx;
  TbxPtr _tbx_idx;
  Kstring _kstr{};
  HtsItrPtr _hts_itr;
};

}  // namespace xoos::io
