#pragma once

#include <xoos/types/fs.h>

namespace xoos::svc {

/**
 * Synopsis:
 * This file contains utility functions for handling files.
 */

void ValidateNonEmptyFile(const fs::path& path);
void ValidateBamAndIndex(const fs::path& path);
void ValidateVcfAndIndex(const fs::path& path);
void ValidateFastaAndIndex(const fs::path& path);
void ValidateBed(const fs::path& path);

/**
 * Given an input bam gets the BAM index file path, either in the form of <input.bam>.bai or <input>.bai. If neither
 * file is present an exception is thrown
 * @param bam_path A sorted and indexed bam file
 * @return a filepath to a bam index file
 */
fs::path GetBamIndexPath(const fs::path& bam_path);

}  // namespace xoos::svc
