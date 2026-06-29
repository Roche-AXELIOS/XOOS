#pragma once

namespace xoos::svc {

/**
 * @brief YC decode methods supported by the Small Variant Caller.
 * @details Supported methods:
 * 1. "none": do not find and decode YC tags
 * 2. "consensus": decode YC tags to infer base types (concordant, simplex, discordant) for consensus reads
 * 3. "split": decode YC tags to split duplex reads into two simplex reads (R1, R2)
 * @note This is used to control how YC tags are handled during feature extraction.
 */
enum class YcDecodeMethod {
  kNone,
  kConsensus,
  kSplit
};

}  // namespace xoos::svc
