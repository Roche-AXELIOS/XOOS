#pragma once
#include <htslib/sam.h>

// NOLINTBEGIN

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Concatenate multiple BAM files into a single BAM file.
 *
 * @param firstfile The first input file opened with hts_open.
 * @param nfn The number of input files.
 * @param fn The array of input file names.
 * @param h The SAM header for the concatenated BAM, can be nullptr.
 * @param outbam The output BAM file name.
 * @param arg_list Additional arguments, currently unused and should be nullptr.
 * @param no_pg Unused, should be set to 1.
 * @return 0 on success, -1 on failure.
 */
int bam_cat(
    samFile* firstfile, int nfn, const char* const* fn, sam_hdr_t* h, const char* outbam, char* arg_list, int no_pg);

#ifdef __cplusplus
}
#endif

// NOLINTEND
