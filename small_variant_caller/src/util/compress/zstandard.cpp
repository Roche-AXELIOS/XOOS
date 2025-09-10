#include "zstandard.h"

#include <zstd.h>

#include <filesystem>
#include <fstream>
#include <iostream>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/vec.h>

#include "util/compress/ios-util.h"

namespace xoos::svc {

void CompressZstd(const fs::path& input_path, const fs::path& output_path, int compression_level, int threads) {
  auto in = OpenInput(input_path, std::ios::binary);
  auto out = OpenOutput(output_path, std::ios::binary);

  const size_t in_buffer_size = ZSTD_CStreamInSize();
  vec<char> in_buffer(in_buffer_size);

  const size_t out_buffer_size = ZSTD_CStreamOutSize();
  vec<char> out_buffer(out_buffer_size);

  auto cctx = ZstdCreateCCtx();

  ZstdSetParameter(cctx.get(), ZSTD_c_compressionLevel, compression_level);
  ZstdSetParameter(cctx.get(), ZSTD_c_checksumFlag, 1);

  if (threads > 1) {
    ZstdSetParameter(cctx.get(), ZSTD_c_nbWorkers, threads);
  }

  const auto to_read = in_buffer_size;
  while (true) {
    // Read compressed data into the input buffer
    const size_t bytes_read = Read(in, in_buffer);
    const auto last_chunk = bytes_read < to_read;
    const auto mode = last_chunk ? ZSTD_e_end : ZSTD_e_continue;

    // Compress the data in chunks until done
    ZSTD_inBuffer input = {in_buffer.data(), bytes_read, 0};
    while (true) {
      ZSTD_outBuffer output = {out_buffer.data(), out_buffer_size, 0};
      const auto remaining = ZstdCompressStream2(cctx.get(), &output, &input, mode);
      out.write(out_buffer.data(), static_cast<std::streamsize>(output.pos));

      if (last_chunk ? remaining == 0 : input.pos == input.size) {
        break;
      }
    }

    if (last_chunk) {
      break;
    }
  }
  RequireIoNotBad(input_path, in);
  RequireIoNotBad(output_path, out);
}

void CompressZstd(const std::string& input_content, const fs::path& output_path, int compression_level) {
  auto out = OpenOutput(output_path, std::ios::binary);

  // Determine the required output buffer size to compress the entire input content, then allocate buffer and compress
  vec<char> buffer_compressed(ZstdCompressBound(input_content.size()));
  const auto compressed_size = ZstdCompress(buffer_compressed, input_content, compression_level);

  out.write(buffer_compressed.data(), static_cast<std::streamsize>(compressed_size));
  RequireIoNotBad(output_path, out);
}

void DecompressZstd(const fs::path& input_path, const fs::path& output_path) {
  auto in = OpenInput(input_path, std::ios::binary);
  auto out = OpenOutput(output_path, std::ios::binary);

  const size_t in_buffer_size = ZSTD_DStreamInSize();
  vec<char> in_buffer(in_buffer_size);

  const size_t out_buffer_size = ZSTD_DStreamOutSize();
  vec<char> out_buffer(out_buffer_size);

  auto dctx = ZstdCreateDCtx();

  size_t last_ret = 0;
  while (true) {
    // Read data into input buffer
    const size_t bytes_read = Read(in, in_buffer);
    if (bytes_read == 0) {
      break;
    }

    // Decompress input buffer in chunks, writing out decompressed data until done
    ZSTD_inBuffer input = {in_buffer.data(), bytes_read, 0};
    while (input.pos < input.size) {
      ZSTD_outBuffer output = {out_buffer.data(), out_buffer_size, 0};
      const auto ret = ZstdDecompressStream(dctx.get(), &output, &input);
      out.write(out_buffer.data(), static_cast<std::streamsize>(output.pos));
      last_ret = ret;
    }
  }

  if (last_ret != 0) {
    throw error::Error("EOF before end of stream: {}", ZSTD_getErrorName(last_ret));
  }

  RequireIoNotBad(input_path, in);
  RequireIoNotBad(output_path, out);
}

std::string DecompressZstd(const fs::path& input_path) {
  auto in = OpenInput(input_path, std::ios::binary);

  const size_t in_buffer_size = ZSTD_DStreamInSize();
  vec<char> in_buffer(in_buffer_size);

  const size_t out_buffer_size = ZSTD_DStreamOutSize();
  vec<char> out_buffer(out_buffer_size);

  auto dctx = ZstdCreateDCtx();

  std::string out;
  out.reserve(ZstdCompressBound(fs::file_size(input_path)));

  size_t last_ret = 0;
  while (true) {
    // Read compressed data into the input buffer
    const size_t bytes_read = Read(in, in_buffer);
    if (bytes_read == 0) {
      break;
    }

    // Decompress input buffer in chunks, writing out decompressed data until done
    ZSTD_inBuffer input = {in_buffer.data(), bytes_read, 0};
    while (input.pos < input.size) {
      ZSTD_outBuffer output = {out_buffer.data(), out_buffer_size, 0};
      const auto ret = ZstdDecompressStream(dctx.get(), &output, &input);
      out.append(out_buffer.data(), output.pos);
      last_ret = ret;
    }
  }

  if (last_ret != 0) {
    throw error::Error("EOF before end of stream: {}", ZSTD_getErrorName(last_ret));
  }

  RequireIoNotBad(input_path, in);
  return out;
}

void ZstdCCtxDeleter::operator()(ZSTD_CCtx* cctx) const {
  if (cctx != nullptr) {
    ZSTD_freeCCtx(cctx);
  }
}

ZstdCCtxPtr ZstdCreateCCtx() {
  ZstdCCtxPtr cctx(ZSTD_createCCtx());
  if (cctx == nullptr) {
    throw error::Error("Failed to create ZSTD_CCtx");
  }
  return cctx;
}

static size_t RequireNoZstdError(const std::string& function_name, const size_t ret) {
  if (ZSTD_isError(ret) == 1) {
    throw error::Error("{} failed: {}", function_name, ZSTD_getErrorName(ret));
  }
  return ret;
}

size_t ZstdSetParameter(ZSTD_CCtx* cctx, ZSTD_cParameter param, int value) {
  return RequireNoZstdError("ZSTD_CCtx_setParameter", ZSTD_CCtx_setParameter(cctx, param, value));
}

size_t ZstdCompressStream2(ZSTD_CCtx* cctx, ZSTD_outBuffer* output, ZSTD_inBuffer* input, ZSTD_EndDirective mode) {
  return RequireNoZstdError("ZSTD_compressStream2", ZSTD_compressStream2(cctx, output, input, mode));
}

void ZstdDCtxDeleter::operator()(ZSTD_DCtx* dctx) const {
  if (dctx != nullptr) {
    ZSTD_freeDCtx(dctx);
  }
}

// Create decompression context with error checking
ZstdDCtxPtr ZstdCreateDCtx() {
  ZstdDCtxPtr dctx(ZSTD_createDCtx());
  if (dctx == nullptr) {
    throw error::Error("Failed to create ZSTD_DCtx");
  }
  return dctx;
}

size_t ZstdDecompressStream(ZSTD_DCtx* dctx, ZSTD_outBuffer* output, ZSTD_inBuffer* input) {
  return RequireNoZstdError("ZSTD_decompressStream", ZSTD_decompressStream(dctx, output, input));
}

size_t ZstdCompressBound(size_t src_size) {
  return RequireNoZstdError("ZSTD_compressBound", ZSTD_compressBound(src_size));
}

size_t ZstdCompress(void* dst, size_t dst_capacity, const void* src, size_t src_size, int compression_level) {
  return RequireNoZstdError("ZSTD_compress", ZSTD_compress(dst, dst_capacity, src, src_size, compression_level));
}

size_t ZstdCompress(vec<char>& dst, const vec<char>& src, int compression_level) {
  return ZstdCompress(dst.data(), dst.size(), src.data(), src.size(), compression_level);
}

size_t ZstdCompress(vec<char>& dst, const std::string& src, int compression_level) {
  return ZstdCompress(dst.data(), dst.size(), src.data(), src.size(), compression_level);
}

}  // namespace xoos::svc
