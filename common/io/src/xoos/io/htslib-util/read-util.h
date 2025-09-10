#pragma once
#include <string>

#include <xoos/types/int.h>

namespace xoos::io {

// Gets the read length by removing the soft clipped bases from the start and end of the read for a given cigar
u32 ApproximateGenomicInsertLength(u32 l_seq, const u32* cigar, u32 n_cigar_op);

// Gets the base for a given string using the htslib direct access via a raw htslib pointer to sequence.
char GetBase(const u8* seq, u32 qpos);

// Convenience function to get the sequence from a read from a raw htslib pointer to sequence.
std::string GetSequence(const u8* seq, u32 start, u32 len);

}  // namespace xoos::io
