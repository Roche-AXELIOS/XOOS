#pragma once

#include <xoos/types/fs.h>

/**
 * @file compress.h
 * @brief Compression and decompression utilities for file operations.
 *
 * This module provides functions for compressing and decompressing files,
 * with automatic format detection based on file extensions.
 */

namespace xoos::svc {

/**
 * @brief Determines if a file is compressed based on its file extension.
 *
 * Checks the file extension to determine if the file is in a compressed format.
 * Supported formats are .gz and .zstd
 *
 * @param path The file path to check
 * @return true if the file appears to be compressed, false otherwise
 *
 * @note This function only checks the file extension, not the actual file content.
 *       It does not verify that the file exists or is actually compressed.
 */
bool IsCompressed(const fs::path& path);

/**
 * @brief Compresses a file from input path to output path.
 *
 * Reads the input file and writes a compressed version to the output path.
 * The compression format is determined by the output file extension.
 *
 * @param input_path Path to the input file to compress
 * @param output_path Path where the compressed file will be written
 *
 * @throws std::runtime_error if input file cannot be read
 * @throws std::runtime_error if output file cannot be written
 * @throws std::runtime_error if compression format is unsupported
 *
 * @note The output directory must exist before calling this function.
 *       The function will overwrite existing files at output_path.
 */
void Compress(const fs::path& input_path, const fs::path& output_path);

/**
 * @brief Compresses string content and writes to output path.
 *
 * Takes string content directly and compresses it to the specified output file.
 * The compression format is determined by the output file extension.
 *
 * @param input_content The string content to compress
 * @param output_path Path where the compressed file will be written
 *
 * @throws std::runtime_error if output file cannot be written
 * @throws std::runtime_error if compression format is unsupported
 *
 * @note Useful for compressing generated content without writing to intermediate files.
 *       The output directory must exist before calling this function.
 */
void Compress(const std::string& input_content, const fs::path& output_path);

/**
 * @brief Decompresses a file from input path to output path.
 *
 * Reads a compressed input file and writes the decompressed content to the output path.
 * The compression format is auto-detected from the input file extension or content.
 *
 * @param input_path Path to the compressed input file
 * @param output_path Path where the decompressed file will be written
 *
 * @throws std::runtime_error if input file cannot be read or is not compressed
 * @throws std::runtime_error if output file cannot be written
 * @throws std::runtime_error if decompression fails due to corrupted data
 *
 * @note The output directory must exist before calling this function.
 *       The function will overwrite existing files at output_path.
 */
void Decompress(const fs::path& input_path, const fs::path& output_path);

/**
 * @brief Decompresses a file and returns the content as a string.
 *
 * Reads a compressed input file and returns the decompressed content as a string.
 * The compression format is auto-detected from the input file extension or content.
 *
 * @param input_path Path to the compressed input file
 * @return The decompressed file content as a string
 *
 * @throws std::runtime_error if input file cannot be read or is not compressed
 * @throws std::runtime_error if decompression fails due to corrupted data
 * @throws std::bad_alloc if the decompressed content is too large for memory
 *
 * @warning Use with caution for large files as the entire content is loaded into memory.
 *          Consider using the file-to-file version for large compressed files.
 */
std::string Decompress(const fs::path& input_path);

}  // namespace xoos::svc
