#pragma once
#include <optional>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>

#include "observations/observations.h"

namespace xoos::cnc {
/**
 * @brief Write BAF values to a file
 * @param baf_values BAF values to write
 * @param out_fname Output file name
 */
void WriteBafValues(const Observations& ref_ads,
                    const Observations& alt_ads,
                    const fs::path& out_fname,
                    const io::CommandLineInfo& command_line_info);

/**
 * @brief Write BAF values to a BigWig file
 * @param ref_ads Reference allele depth observations
 * @param alt_ads Alternate allele depth observations
 * @param out_fname Output BigWig file name
 * @param fai_filename fs::path to FASTA index file containing chromosome lengths
 */
void WriteBafValuesToBigWig(const Observations& ref_ads,
                            const Observations& alt_ads,
                            const fs::path& out_fname,
                            const fs::path& fai_filename);

void WriteBafFiles(const Observations& ref_obvs,
                   const Observations& alt_obvs,
                   const std::optional<fs::path>& baf_out,
                   const std::optional<fs::path>& baf_bw_out,
                   const std::optional<fs::path>& reference_genome_fai_fname,
                   const io::CommandLineInfo& command_line_info);

}  // namespace xoos::cnc
