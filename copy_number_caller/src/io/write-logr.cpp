#include "io/write-logr.h"

#include <fstream>
#include <optional>

#include <xoos/error/error.h>
#include <xoos/io/metadata-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/fs.h>

#include "io/column-names.h"
#include "observations/observations.h"

namespace xoos::cnc {

void WriteLogRFiles(const Observations& logrs,
                    const fs::path& logrs_out,
                    const std::optional<fs::path>& logrs_bw_out,
                    const std::optional<fs::path>& reference_genome_fai_fname,
                    const io::CommandLineInfo& command_line_info) {
  Logging::Info("Writing logrs BED to {}", logrs_out.string());
  std::ofstream ofs(logrs_out);
  if (!ofs.is_open()) {
    throw error::Error("Failed to open logR output file: {}", logrs_out.string());
  }
  logrs.Write(ofs, kColumnLogRatio, false, command_line_info);

  if (logrs_bw_out.has_value() && reference_genome_fai_fname.has_value()) {
    Logging::Info("Writing logrs BigWig to {}", logrs_bw_out.value().string());
    logrs.WriteBigWig(logrs_bw_out.value(), reference_genome_fai_fname.value());
  }
}

void WriteDenoisedLogRFiles(const Observations& logrs,
                            const fs::path& denoised_logrs_out,
                            const io::CommandLineInfo& command_line_info) {
  Logging::Info("Writing logrs BED to {}", denoised_logrs_out.string());
  std::ofstream ofs(denoised_logrs_out);
  if (!ofs.is_open()) {
    throw error::Error("Failed to open denoised logR output file: {}", denoised_logrs_out.string());
  }
  logrs.Write(ofs, kColumnLogRatio, false, command_line_info);
}

}  // namespace xoos::cnc
