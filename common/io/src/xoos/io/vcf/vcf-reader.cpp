#include "xoos/io/vcf/vcf-reader.h"

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-bcf-util.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/util/adopt.h>

namespace xoos::io {

VcfReader::VcfReader(const fs::path& input_vcf_name) {
  // HtsOpen throws if the file cannot be opened; adopt the unique_ptr into the shared handle.
  _input_file = AdoptShared(HtsOpen(input_vcf_name, "r"));

  // Must look for the BGZF EOF marker BEFORE loading the index.
  if (!BcfCheckBgzfEof(_input_file.get())) {
    throw error::Error("Missing BGZF EOF marker; VCF/BCF file may be truncated: {}", input_vcf_name.string());
  }

  const auto compression_type = _input_file->format.compression;
  if (_input_file->format.format == vcf) {
    if (compression_type == bgzf) {
      _tbx_idx = TbxIndexLoad(input_vcf_name);
    }
  } else if (_input_file->format.format == bcf) {
    _hts_idx = BcfIndexLoad(input_vcf_name);
  } else {
    throw error::Error("Not a valid VCF/BCF file type: {}", input_vcf_name.string());
  }

  _hdr = VcfHeader::Read(_input_file);
}

VcfRecordPtr VcfReader::GetNextRecord(const s32 max_unpack) const {
  return VcfRecord::ReadFromFile(_hdr, _input_file, max_unpack);
}

VcfRecordPtr VcfReader::GetNextRecord() const {
  return VcfRecord::ReadFromFile(_hdr, _input_file);
}

std::vector<VcfRecordPtr> VcfReader::GetAllRecords() const {
  std::vector<VcfRecordPtr> records;
  while (true) {
    const auto record = GetNextRecord();
    if (record == nullptr) {
      break;
    }
    records.push_back(record);
  }
  return records;
}

/**
 * @brief Check whether the VCF file opened has an index.
 * @return The index is found or not
 */
bool VcfReader::HasIndex() {
  return _hts_idx != nullptr || _tbx_idx != nullptr;
}

/**
 * @brief Extract contig indexes from the VCF index. If VCF is not indexed, then extract contig indexes from the header.
 * @return Map of contig name to index
 */
StrMap<s32> VcfReader::GetContigIndexes() const {
  auto indexes = _hdr->GetContigIndexes();
  if (_tbx_idx) {
    StrMap<s32> tbx_indexes;
    for (const auto& [chrom, cid] : indexes) {
      const s32 tid = TbxName2Id(_tbx_idx.get(), chrom);
      if (tid >= 0) {
        tbx_indexes[chrom] = tid;
      }
    }
    if (tbx_indexes.empty()) {
      throw error::Error("Empty contig indexes; VCF/BCF index file may be empty or invalid");
    }
    return tbx_indexes;
  }
  return indexes;
}

/**
 * @brief Set the target region for indexed VCF.
 * @param chrom Chromosome
 * @param start Start position, 0-based
 * @param end End position, end-exclusive
 * @return Whether the iterator point to any record
 */
bool VcfReader::SetRegion(const std::string& chrom, const u64 start, const u64 end) {
  if (_input_file->format.format == vcf) {
    if (_input_file->format.compression == bgzf) {
      if (_tbx_idx) {
        const s32 tid = TbxName2Id(_tbx_idx.get(), chrom);
        if (tid >= 0) {
          return SetRegion(tid, start, end);
        }
      } else {
        throw error::Error("VCF file does not have a valid index!");
      }
    } else {
      throw error::Error("VCF file is not bgzip-compressed!");
    }
  } else if (_input_file->format.format == bcf) {
    if (_hts_idx) {
      const s32 tid = _hdr->GetContigId(chrom);
      if (tid >= 0) {
        return SetRegion(tid, start, end);
      }
    } else {
      throw error::Error("BCF file does not have a valid index!");
    }
  } else {
    throw error::Error("File is not a valid VCF/BCF file type!");
  }
  return false;
}

/**
 * @brief Set the target region for indexed VCF.
 * @param tid Chromosome contig index. See `VcfReader::GetContigIndexes`
 * @param start Start position, 0-based
 * @param end End position, end-exclusive
 * @return Whether the iterator point to any record
 */
bool VcfReader::SetRegion(const s32 tid, const u64 start, const u64 end) {
  if (tid >= 0) {
    if (_input_file->format.format == vcf) {
      if (_input_file->format.compression == bgzf) {
        if (_tbx_idx) {
          _hts_itr = TbxItrQueryI(_tbx_idx.get(), tid, static_cast<s64>(start), static_cast<s64>(end));
        } else {
          throw error::Error("VCF file does not have a valid index!");
        }
      } else {
        throw error::Error("VCF file is not bgzip-compressed!");
      }
    } else if (_input_file->format.format == bcf) {
      if (_hts_idx) {
        _hts_itr = BcfItrQueryI(_hts_idx.get(), tid, static_cast<s64>(start), static_cast<s64>(end));
      } else {
        throw error::Error("BCF file does not have a valid index!");
      }
    } else {
      throw error::Error("File is not a valid VCF/BCF file type!");
    }
    return _hts_itr.get() != nullptr;
  }
  return false;
}

/**
 * @brief Return the next record overlapping the target region.
 * @return The next record
 */
VcfRecordPtr VcfReader::GetNextRegionRecord(const s32 unpack) {
  return VcfRecord::ReadFromRegion(_hdr, _input_file, _tbx_idx.get(), _hts_itr.get(), _kstr, unpack);
}

}  // namespace xoos::io
