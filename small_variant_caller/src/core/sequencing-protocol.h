#pragma once

namespace xoos::svc {

/**
 * @brief Sequencing protocols supported by the Small Variant Caller.
 * @details Supported protocols:
 * 1. "duplex": all reads are duplex consensus sequences
 *    - YC tag is used to infer base types (concordant, simplex, discordant)
 * 2. "duplex-simplex": reads are either duplex and simplex reads
 *    - Presence of YC tag indicates a duplex read
 *    - Absence of YC tag indicates a simplex read
 * 3. "umi": all reads are UMI consensus sequences
 * @note This is used to set the ReadType (duplex, simplex, UMI consensus) for each read during feature extraction.
 */
enum class SequencingProtocol {
  kDuplex,
  kDuplexSimplex,
  kUmi
};

/**
 * @brief Checks if the given sequencing protocol involves duplex reads.
 * @param protocol The sequencing protocol to check.
 * @return true if the protocol involves duplex reads, false otherwise.
 */
bool IsDuplexProtocol(SequencingProtocol protocol);

}  // namespace xoos::svc
