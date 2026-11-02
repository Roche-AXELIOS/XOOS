#pragma once
#include <optional>

#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>

#include "observations/observations.h"

namespace xoos::cnc {

void WriteLogRFiles(const Observations& logrs,
                    const fs::path& logrs_out,
                    const std::optional<fs::path>& logrs_bw_out,
                    const std::optional<fs::path>& reference_genome_fai_fname,
                    const io::CommandLineInfo& command_line_info);
void WriteDenoisedLogRFiles(const Observations& logrs,
                            const fs::path& denoised_logrs_out,
                            const io::CommandLineInfo& command_line_info);

}  // namespace xoos::cnc
