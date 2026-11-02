#pragma once

#include <xoos/io/vcf/vcf-header.h>
#include <xoos/io/vcf/vcf-record.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "params.h"
#include "variant.h"

namespace xoos::svc::phaser {

/** @brief Parsed VCF data with the original records preserved for output. */
struct VcfData {
  io::VcfHeaderPtr header;
  vec<io::VcfRecordPtr> records;
  vec<Variant> variants;
  s32 tumor_sample_idx = -1;
};

/**
 * @brief Read a VCF and extract variants for phasing.
 *
 * The tumor sample column is inferred from the VCF header metadata (`tumor_sample`), falling
 * back to the last sample when the metadata is absent. `include_filtered` controls whether
 * non-PASS variants are included, while `rescue_filtered` allows them to be rescued when they
 * phase with a PASS variant.
 *
 * @param vcf_path Input VCF path.
 * @param params Filtering and sample-selection options.
 * @return Parsed VCF header, original records, and variants selected for phasing.
 */
VcfData ReadVcf(const fs::path& vcf_path, const Params& params);

/**
 * @brief Write an output VCF with phasing annotations applied to the tumor sample.
 * @param output_path Destination VCF path.
 * @param vcf_data Parsed input VCF data, including original records and header.
 * @param phase_blocks Phase blocks and haplotype assignments to annotate.
 * @param variants Variants corresponding to the phase-block indices.
 * @param params Merge, filtering, and reference options controlling output.
 */
void WriteVcf(const fs::path& output_path,
              const VcfData& vcf_data,
              const vec<PhaseBlock>& phase_blocks,
              const vec<Variant>& variants,
              const Params& params);

}  // namespace xoos::svc::phaser
