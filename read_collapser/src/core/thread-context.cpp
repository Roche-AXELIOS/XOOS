#include "core/thread-context.h"

#include <optional>

#include <xoos/io/htslib-util/htslib-util.h>

#include "io/alignment-io.h"
#include "io/fastq-writer.h"

namespace xoos::read_collapser {

ThreadContext::ThreadContext(const ReadCollapserOptions& options)
    : _reader(io::OpenAlignmentReader(options.bam_input)), _command_line_info(options.command_line_info) {
}

void ThreadContext::OpenBamWriter(const fs::path& path, const bool write_index) {
  if (!_command_line_info.name.empty()) {
    io::SamHdrAddPgLine(_reader.hdr.get(),
                        {.id = _command_line_info.name,
                         .program_name = _command_line_info.name,
                         .version_number = _command_line_info.version,
                         .command_line = _command_line_info.command_line,
                         .previous_program_id = std::nullopt});
  }
  _bam_writer = OpenAlignmentWriter(path, _reader.hdr.get(), write_index);
}

void ThreadContext::OpenFastqWriter(const fs::path& path, const s32 compression_level) {
  _fastq_writer = OpenGzipFile(path, compression_level);
}

}  // namespace xoos::read_collapser
