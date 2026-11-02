#include "io/fastq-writer.h"

#include <zlib.h>

#include <limits>
#include <ranges>
#include <string>
#include <string_view>

#include <fmt/format.h>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

namespace xoos::read_collapser {

namespace {
// Maximum zlib/gzip compression level (0 = no compression, 9 = maximum).
constexpr s32 kMaxGzipCompressionLevel = 9;

// Write an assembled FASTQ record buffer to the gzip stream, throwing on failure.
void FlushBufferToGzip(gzFile_s* const gz, const std::string_view buf) {
  const s32 ret = gzwrite(gz, buf.data(), static_cast<u32>(buf.size()));  // NOSONAR(S5356) — zlib C API requires void*
  if (ret != static_cast<s32>(buf.size())) {
    s32 err_num = 0;
    const char* const error_string = gzerror(gz, &err_num);
    if (err_num == Z_ERRNO) {
      throw error::Error("Failed to write to file: {}", std::string(error_string));
    }
    throw error::Error("gzip error: {}. Failed to write to file.", std::string(error_string));
  }
}

// Append all duplex consensus auxiliary tags described by `consensus_result` to `buffer`. Requires the per-base
// depth/error vectors to be populated (i.e. --include-consensus-aux-tags was enabled); throws otherwise.
//
// `has_forward`/`has_reverse` report whether the cluster actually contained forward/reverse reads. The per-base array
// tags for a strand with no reads (ad/ae for a reverse-only cluster, bd/be for a forward-only cluster) are all-zero and
// carry no information, so they are omitted. The read-level summary tags (aD/bD/cD, aM/bM/cM, aE/bE/cE) are always
// emitted so downstream tooling sees a fixed set of summary fields regardless of cluster strandedness.
void AppendConsensusAuxTags(std::string& buffer,
                            const ConsensusResult& consensus_result,
                            const bool has_forward,
                            const bool has_reverse) {
  if (!consensus_result.forward_per_base_depth.has_value() || !consensus_result.reverse_per_base_depth.has_value() ||
      !consensus_result.forward_per_base_errors.has_value() || !consensus_result.reverse_per_base_errors.has_value()) {
    throw error::Error("AppendConsensusAuxTags requires per-base depth and error vectors to be populated");
  }

  // Per-base single-strand arrays: depth (ad/bd) and errors (ae/be). Only the strands present in the cluster are
  // emitted; the absent strand's arrays would be uniformly zero.
  auto cast_to_s16 = [](u32 value) -> s16 {
    if (value > std::numeric_limits<s16>::max()) {
      return std::numeric_limits<s16>::max();
    }
    return static_cast<s16>(value);
  };
  if (has_forward) {
    auto downcasted_forward_depth = std::views::transform(*consensus_result.forward_per_base_depth, cast_to_s16);
    auto downcasted_forward_errors = std::views::transform(*consensus_result.forward_per_base_errors, cast_to_s16);
    AppendAuxTag(buffer, "ad", downcasted_forward_depth);
    AppendAuxTag(buffer, "ae", downcasted_forward_errors);
  }
  if (has_reverse) {
    auto downcasted_reverse_depth = std::views::transform(*consensus_result.reverse_per_base_depth, cast_to_s16);
    auto downcasted_reverse_errors = std::views::transform(*consensus_result.reverse_per_base_errors, cast_to_s16);
    AppendAuxTag(buffer, "bd", downcasted_reverse_depth);
    AppendAuxTag(buffer, "be", downcasted_reverse_errors);
  }

  // Per-strand and combined summaries: max depth (aD/bD/cD), min depth (aM/bM/cM), error rate (aE/bE/cE).
  const ConsensusSummary tags = consensus_result.ComputeSummary();
  AppendAuxTag(buffer, "aD", tags.forward.max_depth);
  AppendAuxTag(buffer, "bD", tags.reverse.max_depth);
  AppendAuxTag(buffer, "cD", tags.combined.max_depth);
  AppendAuxTag(buffer, "aM", tags.forward.min_depth);
  AppendAuxTag(buffer, "bM", tags.reverse.min_depth);
  AppendAuxTag(buffer, "cM", tags.combined.min_depth);
  AppendAuxTag(buffer, "aE", tags.forward.error_rate);
  AppendAuxTag(buffer, "bE", tags.reverse.error_rate);
  AppendAuxTag(buffer, "cE", tags.combined.error_rate);
}
}  // namespace

void GzipFileDeleter::operator()(gzFile_s* gz) const {
  if (gz != nullptr) {
    gzclose(gz);
  }
}

void GzipWriteFastq(gzFile_s* const gz,
                    const std::string_view read,
                    const vec<u8>& qual,
                    const std::string_view read_name,
                    std::string& buffer) {
  // clear() retains capacity, so the caller-owned buffer stops reallocating after the first few records.
  buffer.clear();
  buffer += '@';
  buffer += read_name;
  buffer += '\n';
  buffer += read;
  buffer += "\n+\n";
  for (const u8 q : qual) {
    buffer += ToPlus33Ascii(q);
  }
  buffer += '\n';
  FlushBufferToGzip(gz, buffer);
}

void GzipWriteFastq(gzFile_s* const gz,
                    const ConsensusResult& consensus_result,
                    const std::string_view read_name,
                    std::string& buffer,
                    const bool has_forward,
                    const bool has_reverse) {
  // clear() retains capacity, so the caller-owned buffer stops reallocating after the first few records.
  buffer.clear();
  buffer += '@';
  buffer += read_name;
  // The four per-base vectors are populated together (or not at all) by MajorityVotingConsensus, so a single
  // presence check gates the whole aux tag block.
  if (consensus_result.forward_per_base_depth.has_value()) {
    AppendConsensusAuxTags(buffer, consensus_result, has_forward, has_reverse);
  }
  buffer += '\n';
  buffer += consensus_result.sequence;
  buffer += "\n+\n";
  for (const u8 q : consensus_result.quality_scores) {
    buffer += ToPlus33Ascii(q);
  }
  buffer += '\n';
  FlushBufferToGzip(gz, buffer);
}

GzipFilePtr OpenGzipFile(const fs::path& path, s32 compression_level) {
  create_directories(path.parent_path());
  if (compression_level > kMaxGzipCompressionLevel) {
    throw std::runtime_error("Invalid compression level: " + std::to_string(compression_level) + ", must be [0, 9]");
  }
  std::string write_mode = compression_level == 0 ? "wT" : "w" + std::to_string(compression_level);
  GzipFilePtr gz(gzopen(path.c_str(), write_mode.c_str()));
  if (!gz) {
    throw error::Error("Failed to open gzip file: '{}'", path.string());
  }
  return gz;
}

}  // namespace xoos::read_collapser
