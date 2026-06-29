#pragma once

#include <xoos/types/fs.h>

namespace xoos::svc {

/**
 * Given an input bam gets the BAM index file path, either in the form of <input.bam>.bai or <input>.bai. If neither
 * file is present an exception is thrown
 * @param bam_path A sorted and indexed bam file
 * @return a filepath to a bam index file
 */
fs::path GetBamIndexPath(const fs::path& bam_path);

/**
 * @brief Create the parent directory for a file path if it does not exist.
 * @param file_path File path to ensure parent directory exists
 */
void CreateParentDirectoryIfNotExists(const fs::path& file_path);

}  // namespace xoos::svc
