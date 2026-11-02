#pragma once

#include <xoos/io/alignment-reader.h>
#include <xoos/io/metadata-util.h>
#include <xoos/types/fs.h>
#include <xoos/types/int.h>

#include "core/read-collapser-options.h"
#include "io/alignment-io.h"
#include "io/fastq-writer.h"
#include "metrics/metrics.h"

namespace xoos::read_collapser {

/**
 * Per-thread resources for the read collapser. Each processing thread owns one `ThreadContext`: its own alignment
 * reader, an optional BAM writer and/or FASTQ writer, and a `Metrics` accumulator. Owning these by value (rather than
 * referencing driver-owned parallel vectors) keeps each thread's resources together and lets the BAM writer flush and
 * index by RAII when the context is destroyed. The markdup driver relies on this: it snapshots each context's metrics,
 * then destroys the contexts (by scope) before Phase 2 re-reads the temp BAMs.
 *
 * Construction opens the reader; the BAM and FASTQ writers are opened on demand via `OpenBamWriter` / `OpenFastqWriter`
 * so each mode requests only the outputs it needs. Callers read from `Reader()`, write via `BamWriter()` /
 * `FastqWriter()`, and accumulate into `GetMetrics()`. Drive it with the free `ProcessRegionInBatches` (see
 * region-processor.h).
 *
 * When a writer has not been opened, `BamWriter().bam` (respectively `FastqWriter()`) is null; callers must gate writes
 * on the relevant option.
 */
class ThreadContext {
 public:
  // Opens the per-thread alignment reader from `options.bam_input`. Retains the command-line info (by value) so a PG
  // line can be added to any BAM writer opened later.
  explicit ThreadContext(const ReadCollapserOptions& options);

  // Opens a BAM writer at @p path, adding the module PG line to the reader's header. When @p write_index is true an
  // index is written alongside. The writer is flushed and indexed by RAII on destruction.
  void OpenBamWriter(const fs::path& path, bool write_index);

  // Opens a gzip FASTQ writer at @p path with the given @p compression_level.
  void OpenFastqWriter(const fs::path& path, s32 compression_level);

  // The per-thread alignment reader.
  const io::AlignmentReader& Reader() const {
    return _reader;
  }

  // The BAM writer; its `bam` handle is null until OpenBamWriter is called.
  const AlignmentWriter& BamWriter() const {
    return _bam_writer;
  }

  // The FASTQ writer; null until OpenFastqWriter is called.
  const GzipFilePtr& FastqWriter() const {
    return _fastq_writer;
  }

  // The per-thread metrics accumulator.
  Metrics& GetMetrics() {
    return _metrics;
  }

  const Metrics& GetMetrics() const {
    return _metrics;
  }

 private:
  // One alignment reader per thread so concurrent threads do not share a file handle.
  io::AlignmentReader _reader;
  // BAM writer; empty (null `bam`) until OpenBamWriter is called. Flushed and indexed by RAII on destruction.
  AlignmentWriter _bam_writer;
  // FASTQ writer; null until OpenFastqWriter is called.
  GzipFilePtr _fastq_writer;
  // Per-thread metrics accumulated during processing and merged by the driver.
  Metrics _metrics;
  // Command-line info for the PG line, kept by value to avoid coupling to the caller's `options` lifetime.
  io::CommandLineInfo _command_line_info;
};

}  // namespace xoos::read_collapser
