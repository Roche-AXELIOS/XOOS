#pragma once

#include <zlib.h>

#include <memory>

#include <xoos/types/fs.h>

/**
 * @file gzip.h
 * @brief Gzip compression and decompression utilities.
 *
 * This module provides high-level functions for compressing and decompressing
 * files using the gzip compression algorithm, along with low-level wrapper
 * functions for the zlib C API.
 */

namespace xoos::svc {

/**
 * @brief Compresses a file using gzip compression.
 *
 * Reads the input file in chunks and writes a compressed version using
 * the specified buffer size and compression level.
 *
 * @param input_path Path to the input file to compress
 * @param output_path Path where the compressed file will be written
 * @param buffer_size Size of the buffer for reading chunks (affects memory usage)
 * @param compression_level Compression level (0-9, where 9 = best compression)
 *
 * @throws std::runtime_error if files cannot be opened or compression fails
 *
 * @note Larger buffer sizes may improve performance but use more memory.
 *       Default strategy (Z_DEFAULT_STRATEGY) is used for compression.
 */
void CompressGzip(const fs::path& input_path, const fs::path& output_path, size_t buffer_size, int compression_level);

/**
 * @brief Compresses string content using gzip compression.
 *
 * Takes string content directly and compresses it to the specified output file.
 * Uses default strategy for compression.
 *
 * @param input_content The string content to compress
 * @param output_path Path where the compressed file will be written
 * @param compression_level Compression level (0-9, where 9 = best compression)
 *
 * @throws std::runtime_error if output file cannot be written or compression fails
 *
 * @note Suitable for smaller content that fits comfortably in memory.
 *       For large content, consider using the file-based version with buffering.
 */
void CompressGzip(const std::string& input_content, const fs::path& output_path, int compression_level);

/**
 * @brief Decompresses a gzip-compressed file.
 *
 * Reads a compressed input file in chunks and writes the decompressed content
 * to the output file using the specified buffer size.
 *
 * @param input_path Path to the compressed input file
 * @param output_path Path where the decompressed file will be written
 * @param buffer_size Size of the buffer for reading chunks (affects memory usage)
 *
 * @throws std::runtime_error if files cannot be opened or decompression fails
 *
 * @note Larger buffer sizes may improve performance but use more memory.
 *       The function handles EOF detection automatically.
 */
void DecompressGzip(const fs::path& input_path, const fs::path& output_path, size_t buffer_size);

/**
 * @brief Decompresses a gzip-compressed file and returns content as string.
 *
 * Reads a compressed input file and returns the decompressed content as a string.
 * Uses estimated compression ratio for memory reservation.
 *
 * @param input_path Path to the compressed input file
 * @param buffer_size Size of the buffer for reading chunks
 * @return The decompressed content as a string
 *
 * @throws std::runtime_error if file cannot be opened or decompression fails
 * @throws std::bad_alloc if decompressed content is too large for memory
 *
 * @warning Loads entire decompressed content into memory.
 *          Use file-to-file version for very large compressed files.
 * @note Reserves memory based on estimated 3:1 compression ratio.
 */
std::string DecompressGzip(const fs::path& input_path, size_t buffer_size);

/**
 * @brief Custom deleter for gzFile handles.
 *
 * Ensures proper cleanup of gzip file handles by calling gzclose.
 */
struct GzFileDeleter {
  void operator()(gzFile file) const;
};

using GzFilePtr = std::unique_ptr<std::remove_pointer_t<gzFile>, GzFileDeleter>;

/**
 * @brief Opens a gzip file with error checking.
 *
 * Opens a gzip file for reading or writing with the specified mode.
 * Returns a smart pointer that automatically closes the file.
 *
 * @param path Path to the gzip file
 * @param mode Opening mode ("rb" for read, "wb" for write, etc.)
 * @return Smart pointer to the opened gzip file
 *
 * @throws std::runtime_error if the file cannot be opened
 *
 * @note The file is automatically closed when the smart pointer is destroyed.
 *       Common modes: "rb" (read binary), "wb" (write binary).
 */
GzFilePtr GzOpen(const fs::path& path, const std::string& mode);

/**
 * @brief Checks for gzip errors and throws if any are found.
 *
 * Examines the error state of a gzip file and throws an exception
 * if any errors are detected.
 *
 * @param file The gzip file to check
 *
 * @throws std::runtime_error if the file has encountered an error
 *
 * @note Should be called after gzip operations to detect errors.
 *       Only throws if the error code is not Z_OK.
 */
void GzError(gzFile file);

/**
 * @brief Writes data to a gzip file with error checking.
 *
 * Writes the specified data to a gzip file and checks for errors.
 *
 * @param file The gzip file to write to
 * @param data Pointer to the data to write
 * @param size Number of bytes to write
 *
 * @throws std::runtime_error if the write operation fails
 *
 * @note The entire data buffer is written in a single operation.
 *       For large data, consider writing in chunks.
 */
void GzWrite(gzFile file, const char* data, size_t size);

/**
 * @brief Sets compression parameters for a gzip file.
 *
 * Configures the compression level and strategy for an open gzip file.
 * Must be called before writing data.
 *
 * @param file The gzip file to configure
 * @param level Compression level (0-9, where 9 = best compression)
 * @param strategy Compression strategy (e.g., Z_DEFAULT_STRATEGY)
 *
 * @throws std::runtime_error if parameter setting fails
 *
 * @note Must be called on a file opened for writing before any data is written.
 *       Common strategies: Z_DEFAULT_STRATEGY, Z_FILTERED, Z_HUFFMAN_ONLY.
 */
void GzSetParams(gzFile file, int level, int strategy);

/**
 * @brief Reads data from a gzip file with error checking.
 *
 * Reads up to the specified number of bytes from a gzip file.
 * Returns the actual number of bytes read.
 *
 * @param file The gzip file to read from
 * @param buffer Buffer to read data into
 * @param size Maximum number of bytes to read
 * @return Number of bytes actually read
 *
 * @throws std::runtime_error if the read operation fails
 *
 * @note Returns fewer bytes than requested when EOF is reached.
 *       The buffer must be pre-allocated to at least 'size' bytes.
 */
int GzRead(gzFile file, char* buffer, size_t size);

}  // namespace xoos::svc
