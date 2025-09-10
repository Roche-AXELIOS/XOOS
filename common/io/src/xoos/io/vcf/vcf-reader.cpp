#include "xoos/io/vcf/vcf-reader.h"

#include <stdexcept>

#include <htslib/bgzf.h>

namespace xoos::io {

void TbxDeleter::operator()(tbx_t* tbx) const {
  tbx_destroy(tbx);
}

VcfReader::VcfReader(const std::string& input_vcf_name) {
  htsFile* fp = bcf_open(input_vcf_name.c_str(), "r");
  if (fp == nullptr) {
    throw std::runtime_error("Could not open VCF/BCF file " + input_vcf_name);
  }
  // NOTE: bcf_close defines hts_close, so this is the proper pointer
  _input_file = HtsFileSharedPtr{fp, hts_close};
  if (_input_file == nullptr) {
    throw std::runtime_error("Could not open VCF/BCF file " + input_vcf_name);
  }

  auto compression_type = _input_file->format.compression;
  if (compression_type == bgzf) {
    BGZF* bgzf = hts_get_bgzfp(_input_file.get());
    if ((bgzf != nullptr) && bgzf_check_EOF(bgzf) == 0) {
      // must look for EOF marker BEFORE loading the index
      throw std::runtime_error("Missing BGZF EOF marker; VCF/BCF file may be truncated: " + input_vcf_name);
    }
  }

  if (_input_file->format.format == vcf) {
    if (compression_type == bgzf) {
      _tbx_idx = TbxPtr{tbx_index_load3(input_vcf_name.c_str(), nullptr, HTS_IDX_SILENT_FAIL)};
    }
  } else if (_input_file->format.format == bcf) {
    _hts_idx = HtsIdxPtr{bcf_index_load3(input_vcf_name.c_str(), nullptr, HTS_IDX_SILENT_FAIL)};
  } else {
    throw std::runtime_error("Not a valid VCF/BCF file type: " + input_vcf_name);
  }

  _hdr = VcfHeader::Read(_input_file);
}

VcfRecordPtr VcfReader::GetNextRecord(int max_unpack) const {
  return VcfRecord::ReadFromFile(_hdr, _input_file, max_unpack);
}

VcfRecordPtr VcfReader::GetNextRecord() const {
  return VcfRecord::ReadFromFile(_hdr, _input_file);
}

std::vector<VcfRecordPtr> VcfReader::GetAllRecords() const {
  std::vector<VcfRecordPtr> records;
  while (true) {
    auto record = GetNextRecord();
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
std::map<std::string, s32> VcfReader::GetContigIndexes() {
  auto indexes = _hdr->GetContigIndexes();
  if (_tbx_idx) {
    std::map<std::string, s32> tbx_indexes;
    for (const auto& [chrom, cid] : indexes) {
      const int tid = tbx_name2id(_tbx_idx.get(), chrom.c_str());
      if (tid >= 0) {
        tbx_indexes[chrom] = tid;
      }
    }
    if (tbx_indexes.empty()) {
      throw std::runtime_error("Empty contig indexes; VCF/BCF index file may be empty or invalid");
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
bool VcfReader::SetRegion(const std::string& chrom, u64 start, u64 end) {
  if (_input_file->format.format == vcf) {
    if (_input_file->format.compression == bgzf) {
      if (_tbx_idx) {
        const int tid = tbx_name2id(_tbx_idx.get(), chrom.c_str());
        if (tid >= 0) {
          return SetRegion(tid, start, end);
        }
      } else {
        throw std::runtime_error("VCF file does not have a valid index!");
      }
    } else {
      throw std::runtime_error("VCF file is not bgzip-compressed!");
    }
  } else if (_input_file->format.format == bcf) {
    if (_hts_idx) {
      const int tid = _hdr->GetContigId(chrom);
      if (tid >= 0) {
        return SetRegion(tid, start, end);
      }
    } else {
      throw std::runtime_error("BCF file does not have a valid index!");
    }
  } else {
    throw std::runtime_error("File is not a valid VCF/BCF file type!");
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
bool VcfReader::SetRegion(int tid, u64 start, u64 end) {
  if (tid >= 0) {
    if (_input_file->format.format == vcf) {
      if (_input_file->format.compression == bgzf) {
        if (_tbx_idx) {
          _hts_itr = HtsItrPtr(tbx_itr_queryi(_tbx_idx.get(), tid, start, end));
        } else {
          throw std::runtime_error("VCF file does not have a valid index!");
        }
      } else {
        throw std::runtime_error("VCF file is not bgzip-compressed!");
      }
    } else if (_input_file->format.format == bcf) {
      if (_hts_idx) {
        _hts_itr = HtsItrPtr(bcf_itr_queryi(_hts_idx.get(), tid, start, end));
      } else {
        throw std::runtime_error("BCF file does not have a valid index!");
      }
    } else {
      throw std::runtime_error("File is not a valid VCF/BCF file type!");
    }
    return _hts_itr.get() != nullptr;
  }
  return false;
}

/**
 * @brief Return the next record overlapping the target region.
 * @return The next record
 */
VcfRecordPtr VcfReader::GetNextRegionRecord(const int unpack) {
  return VcfRecord::ReadFromRegion(_hdr, _input_file, _tbx_idx.get(), _hts_itr.get(), _kstr, unpack);
}

}  // namespace xoos::io
