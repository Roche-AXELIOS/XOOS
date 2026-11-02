#pragma once

#include <htslib/sam.h>

#include <xoos/types/int.h>

// Test fixture provenance (common/io/tests/resources/bam/):
//
// The aligner-detection unit tests use two header-only BAM fixtures (no read records) extracted from real
// Genia-aligned BAMs and trimmed to a representative subset of @PG lines (the aligner record, a few trailing
// samtools records to exercise multi-@PG scanning) plus a couple of @SQ lines:
//
//   genia_bwa_header.bam     aligner @PG: ID:bwa PN:bwa VN:0.7.18 CL:bwa mem ...
//                            source: functional_tests/str_caller/resources/bam_input/
//                                    str_HG002_sbx_fast_bwa_with_chrY.bam
//   genia_giraffe_header.bam aligner @PG: ID:0 PN:vg
//                            source: functional_tests/str_caller/resources/bam_input/
//                                    str_HG003_sbx_d_giraffe_with_chrY.bam
//
// The giraffe fixture's aligner record has a numeric ID (ID:0) and records the program only in PN:vg. This is
// why DetectAlignerFromHeader must read the PN tag and cannot rely on the @PG ID alone.

namespace xoos::io {

// Aligner detected from a BAM file, collapsed into the alignment families XOOS distinguishes downstream.
// The XOOS pipeline aligners (bwa, minimap2, vg_giraffe, pb_fq2bam, pb_minimap2, pb_giraffe) map onto these:
// linear bwa-style alignments (bwa / bwa-mem2 / parabricks fq2bam), graph alignments (vg giraffe / parabricks
// giraffe), and minimap2 alignments (minimap2 / parabricks minimap2).
// kUnknown is returned by DetectAlignerFromHeader when no recognized aligner is found.
enum class Aligner : u8 {
  // Returned when the aligner cannot be determined from the BAM header.
  kUnknown,
  kBwa,
  kGiraffe,
  kMinimap2
};

/**
 * @brief Infer the aligner family from a BAM header's @PG program lines.
 * @details Scans the @PG lines for known aligner program names, checking the program name (PN) tag, then the
 *          program record ID, then the command line (CL). Matching is case-insensitive and substring-based:
 *          `giraffe`/`vg` map to giraffe; `fq2bam`, `bwa` and `bwa-mem2` map to bwa; `minimap2` maps to
 *          minimap2. Parabricks records the subcommand in the PN/ID tags (`pbrun fq2bam`/`pbrun giraffe`), and
 *          these tokens are matched there and, defensively, in the command line. Returns `kUnknown` if no
 *          recognized aligner is found.
 */
Aligner DetectAlignerFromHeader(sam_hdr_t* header);

}  // namespace xoos::io
