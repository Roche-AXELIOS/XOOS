#pragma once

#include <zstd.h>

#include <memory>
#include <string>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

/**
 * @file zstandard.h
 * @brief Zstandard compression and decompression utilities.
 *
 * This module provides high-level functions for compressing and decompressing
 * files using the Zstandard (zstd) compression algorithm, along with low-level
 * wrapper functions for the zstd C API.
 */

namespace xoos::svc {

/**
 * @brief Custom deleter for ZSTD_CCtx compression context.
 */
struct ZstdCCtxDeleter {
  void operator()(ZSTD_CCtx* cctx) const;
};

/**
 * @brief Custom deleter for ZSTD_DCtx decompression context.
 */
struct ZstdDCtxDeleter {
  void operator()(ZSTD_DCtx* dctx) const;
};

using ZstdCCtxPtr = std::unique_ptr<ZSTD_CCtx, ZstdCCtxDeleter>;
using ZstdDCtxPtr = std::unique_ptr<ZSTD_DCtx, ZstdDCtxDeleter>;

/**
 * @brief Compresses a file using Zstandard compression.
 *
 * Reads the input file and writes a compressed version using streaming compression.
 * Supports multi-threading for improved performance on large files.
 *
 * @param input_path Path to the input file to compress
 * @param output_path Path where the compressed file will be written
 * @param compression_level Compression level (1-22, higher = better compression)
 * @param threads Number of worker threads to use (1 = single-threaded)
 *
 * @throws std::runtime_error if files cannot be opened or compression fails
 *
 * @note Uses streaming compression for memory efficiency with large files.
 *       Includes checksum verification for data integrity.
 */
void CompressZstd(const fs::path& input_path, const fs::path& output_path, int compression_level, int threads);

/**
 * @brief Compresses string content using Zstandard compression.
 *
 * Takes string content directly and compresses it to the specified output file.
 * Uses simple compression (not streaming) for small to medium-sized content.
 *
 * @param input_content The string content to compress
 * @param output_path Path where the compressed file will be written
 * @param compression_level Compression level (1-22, higher = better compression)
 *
 * @throws std::runtime_error if output file cannot be written or compression fails
 *
 * @note Loads entire content into memory - suitable for smaller data.
 *       For large content, consider using the file-based version.
 */
void CompressZstd(const std::string& input_content, const fs::path& output_path, int compression_level);

/**
 * @brief Decompresses a Zstandard-compressed file.
 *
 * Reads a compressed input file and writes the decompressed content to the output file.
 * Uses streaming decompression for memory efficiency.
 *
 * @param input_path Path to the compressed input file
 * @param output_path Path where the decompressed file will be written
 *
 * @throws std::runtime_error if files cannot be opened, decompression fails, or data is corrupted
 *
 * @note Validates that the entire stream was properly decompressed.
 *       Uses streaming decompression for memory efficiency with large files.
 */
void DecompressZstd(const fs::path& input_path, const fs::path& output_path);

/**
 * @brief Decompresses a Zstandard-compressed file and returns content as string.
 *
 * Reads a compressed input file and returns the decompressed content as a string.
 * Uses streaming decompression but accumulates result in memory.
 *
 * @param input_path Path to the compressed input file
 * @return The decompressed content as a string
 *
 * @throws std::runtime_error if file cannot be opened, decompression fails, or data is corrupted
 * @throws std::bad_alloc if decompressed content is too large for memory
 *
 * @warning Loads entire decompressed content into memory.
 *          Use file-to-file version for very large compressed files.
 */
std::string DecompressZstd(const fs::path& input_path);

// Low-level zstd wrapper functions

/**
 * @brief Creates a Zstandard compression context with error checking.
 */
ZstdCCtxPtr ZstdCreateCCtx();

/**
 * @brief Sets a compression parameter on the context.
 */
size_t ZstdSetParameter(ZSTD_CCtx* cctx, ZSTD_cParameter param, int value);

/**
 * @brief Performs streaming compression.
 */
size_t ZstdCompressStream2(ZSTD_CCtx* cctx, ZSTD_outBuffer* output, ZSTD_inBuffer* input, ZSTD_EndDirective mode);

/**
 * @brief Creates a Zstandard decompression context with error checking.
 */
ZstdDCtxPtr ZstdCreateDCtx();

/**
 * @brief Performs streaming decompression.
 */
size_t ZstdDecompressStream(ZSTD_DCtx* dctx, ZSTD_outBuffer* output, ZSTD_inBuffer* input);

/**
 * @brief Calculates maximum compressed size for given input size.
 */
size_t ZstdCompressBound(size_t src_size);

/**
 * @brief Simple compression function (non-streaming).
 */
size_t ZstdCompress(void* dst, size_t dst_capacity, const void* src, size_t src_size, int compression_level);

/**
 * @brief Simple compression function for vectors.
 */
size_t ZstdCompress(vec<char>& dst, const vec<char>& src, int compression_level);

/**
 * @brief Simple compression function for strings.
 */
size_t ZstdCompress(vec<char>& dst, const std::string& src, int compression_level);

}  // namespace xoos::svc
