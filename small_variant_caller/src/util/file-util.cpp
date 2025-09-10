#include "file-util.h"

#include <zlib.h>
#include <zstd.h>

#include <filesystem>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <xoos/error/error.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

namespace xoos::svc {

/**
 * @brief Validates that a given file exists and is not empty.
 * @param path Path to the file to validate.
 * @throws std::runtime_error If the file does not exist or is empty.
 */
void ValidateNonEmptyFile(const fs::path& path) {
  if (fs::file_size(path) == 0) {
    throw error::Error("File is empty: {}", path);
  }
}

/**
 * @brief Validates that a given VCF file is readable, correctly formatted, and has a valid index.
 * @param path Path to the VCF file to validate.
 * @throws std::runtime_error Any issue encountered during validation
 */
void ValidateVcfAndIndex(const fs::path& path) {
  io::HtsFilePtr vcf_file(hts_open(path.c_str(), "r"));
  if (vcf_file == nullptr) {
    throw error::Error("Failed to open VCF file: '{}'", path);
  }
  if (vcf_file->format.format != vcf) {
    throw error::Error("Input file is not in VCF format: '{}'", path);
  }
  {
    BGZF* bgzf = hts_get_bgzfp(vcf_file.get());
    if (bgzf == nullptr) {
      throw error::Error("VCF file is not bgzip-compressed: {}", path);
    }
    if (bgzf_check_EOF(bgzf) == 0) {
      throw error::Error("VCF file may be truncated: {}", path);
    }
  }
  const auto hdr = io::BcfHeaderPtr(bcf_hdr_read(vcf_file.get()), bcf_hdr_destroy);
  if (hdr == nullptr) {
    throw error::Error("Failed to read VCF header: {}", path);
  }
  auto record = io::BcfRecordPtr{bcf_init(), bcf_destroy};
  if (record == nullptr) {
    throw error::Error("Failed to allocate memory for VCF record");
  }
  const int ret = bcf_read(vcf_file.get(), hdr.get(), record.get());
  if (ret < 0) {
    throw error::Error("VCF file is empty: '{}'", path);
  }
  const io::HtsIdxPtr idx(hts_idx_load(path.c_str(), HTS_FMT_TBI));
  if (idx == nullptr) {
    throw error::Error("Failed to open VCF index file: '{}'", path);
  }
}

/**
 * @brief Validates that a given BAM file is readable, correctly formatted, and has a valid index.
 * @param path Path to the BAM file to validate.
 * @throws std::runtime_error Any issue encountered during validation
 */
void ValidateBamAndIndex(const fs::path& path) {
  io::HtsFilePtr bam_file(hts_open(path.c_str(), "r"));
  if (bam_file == nullptr) {
    throw error::Error("Failed to open BAM file: '{}'", path);
  }
  if (bam_file->format.format != bam) {
    throw error::Error("Input file is not in BAM format: '{}'", path);
  }
  {
    BGZF* bgzf = hts_get_bgzfp(bam_file.get());
    if (bgzf == nullptr) {
      throw error::Error("BAM file is not bgzip-compressed: {}", path);
    }
    if (bgzf_check_EOF(bgzf) == 0) {
      throw error::Error("BAM file may be truncated: {}", path);
    }
  }
  const io::SamHdrPtr hdr(sam_hdr_read(bam_file.get()));
  if (hdr == nullptr) {
    throw error::Error("Failed to read BAM header: {}", path);
  }
  auto record = io::Bam1Ptr(bam_init1());
  if (record == nullptr) {
    throw error::Error("Failed to allocate memory for BAM record");
  }
  const int ret = sam_read1(bam_file.get(), hdr.get(), record.get());
  if (ret == -1) {
    // EOF reached
    throw error::Error("BAM file is empty: '{}'", path);
  }
  if (ret < -1) {
    // Error reading the BAM file
    throw error::Error("Error reading BAM file: '{}'", path);
  }
  const io::HtsIdxPtr idx(sam_index_load(bam_file.get(), path.c_str()));
  if (idx == nullptr) {
    throw error::Error("Failed to open BAM index file");
  }
}

/**
 * @brief Validates that a given FASTA file is readable, correctly formatted, and has a valid index.
 * @param path Path to the FASTA file to validate.
 * @throws std::runtime_error Any issue encountered during validation
 */
void ValidateFastaAndIndex(const fs::path& path) {
  io::HtsFilePtr fasta_file(hts_open(path.c_str(), "r"));
  if (fasta_file == nullptr) {
    throw error::Error("Failed to open FASTA file: '{}'", path);
  }
  if (fasta_file->format.format != fasta_format) {
    throw error::Error("Input file is not in FASTA format: '{}'", path);
  }
  const io::FaIdxPtr fai(fai_load(path.c_str()));
  if (fai == nullptr) {
    throw error::Error("Failed to open FASTA index file: '{}'", path);
  }
}

/**
 * @brief Validates that a given BED file is readable.
 * @param path Path to the BED file to validate.
 * @throws std::runtime_error Any issue encountered during validation
 */
void ValidateBed(const fs::path& path) {
  io::HtsFilePtr bed_file(hts_open(path.c_str(), "r"));
  if (bed_file == nullptr) {
    throw error::Error("Failed to open BED file: '{}'", path);
  }
  if (bed_file->format.format != bed) {
    throw error::Error("Input file is not in BED format: '{}'", path);
  }
}

fs::path GetBamIndexPath(const fs::path& bam_path) {
  // NOTE: CRAM support for feature extraction not yet implemented
  const std::string index_extension = (bam_path.extension() == ".cram") ? ".crai" : ".bai";
  auto index_path = fs::path{bam_path.string() + index_extension};
  if (fs::exists(index_path)) {
    return index_path;
  }
  auto index_fallback_path = fs::path{bam_path}.replace_extension(index_extension);
  if (fs::exists(index_fallback_path)) {
    return index_fallback_path;
  }
  throw error::Error("No {} index exists for input {}", index_extension, bam_path);
}

}  // namespace xoos::svc
