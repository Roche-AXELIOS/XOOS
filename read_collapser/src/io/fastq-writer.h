#pragma once

#include <zlib.h>

#include <concepts>
#include <ranges>
#include <string>
#include <string_view>

#include <xoos/types/fs.h>
#include <xoos/types/int.h>
#include <xoos/types/vec.h>

#include "consensus/consensus-result.h"

namespace xoos::read_collapser {

// Maximum Phred quality value representable in the FASTQ Phred+33 encoding (ASCII '~').
constexpr u8 kMaxPhredPlus33 = 93;
// ASCII offset applied to Phred quality values in the FASTQ Phred+33 encoding (ASCII '!').
constexpr u8 kPhred33Offset = 33;

// Cap at kMaxPhredPlus33 and add kPhred33Offset to fit in ASCII printable range
template <std::unsigned_integral T>
static inline char ToPlus33Ascii(const T val) {
  return static_cast<char>(std::min(val, static_cast<T>(kMaxPhredPlus33)) + static_cast<T>(kPhred33Offset));
}

// A single-value BAM tag is a string (SAM Z), or a scalar integer/float (SAM i/f).
template <typename T>
concept ConvertibleToBamSingleValueTag =
    std::convertible_to<std::remove_cvref_t<T>, std::string_view> || std::integral<std::remove_cvref_t<T>> ||
    std::floating_point<std::remove_cvref_t<T>>;

// An array BAM tag is a range of integers or floats (SAM B). Strings are excluded so they render as SAM Z rather than a
// comma-separated B array.
template <typename T>
concept ConvertibleToBamArrayTag =
    !std::convertible_to<std::remove_cvref_t<T>, std::string_view> && std::ranges::range<std::remove_cvref_t<T>> &&
    (std::integral<std::ranges::range_value_t<std::remove_cvref_t<T>>> ||
     std::floating_point<std::ranges::range_value_t<std::remove_cvref_t<T>>>);

template <typename T>
concept ConvertibleToBamTag = ConvertibleToBamSingleValueTag<T> || ConvertibleToBamArrayTag<T>;

// Return the single-character SAM auxiliary field type code for tag value type `T` (Z string, A char, i integer,
// f float, or B array).
template <typename T>
consteval char ToBamType() {
  using U = std::remove_cvref_t<T>;

  if constexpr (std::convertible_to<U, std::string_view>) {
    return 'Z';
  } else if constexpr (std::same_as<U, char>) {
    return 'A';
  } else if constexpr (std::integral<U>) {
    return 'i';
  } else if constexpr (std::floating_point<U>) {
    return 'f';
  } else if constexpr (std::ranges::range<U>) {
    return 'B';
  } else {
    static_assert(false, "Unsupported type for BAM tag");
  }
}

// Return the SAM B-array element type code for element type `T` (e.g. 's' for s16, 'S' for u16, 'f' for f32).
template <typename T>
consteval char ToBamArrayElementType() {
  using U = std::remove_cvref_t<T>;

  if constexpr (std::same_as<U, s8>) {
    return 'c';
  } else if constexpr (std::same_as<U, u8>) {
    return 'C';
  } else if constexpr (std::same_as<U, s16>) {
    return 's';
  } else if constexpr (std::same_as<U, u16>) {
    return 'S';
  } else if constexpr (std::same_as<U, s32>) {
    return 'i';
  } else if constexpr (std::same_as<U, u32>) {
    return 'I';
  } else if constexpr (std::same_as<U, f32>) {
    return 'f';
  } else {
    static_assert(false, "Unsupported type for BAM array element");
  }
}

// Append a single SAM auxiliary field to `buffer`, formatted as "\t<tag>:<type>:<value>". Strings render as SAM Z,
// scalars as SAM i/f, and ranges as SAM B arrays with comma-separated elements.
template <ConvertibleToBamTag T>
void AppendAuxTag(std::string& buffer, const std::string_view tag, const T& value) {
  buffer += '\t';
  buffer += tag;
  buffer += ':';

  using U = std::remove_cvref_t<T>;

  if constexpr (std::convertible_to<U, std::string_view>) {
    // A string renders as a single SAM Z field. This is checked before the array branch because a string is also a
    // range<char>, which would otherwise be encoded as a comma-separated B array.
    buffer += ToBamType<U>();
    buffer += ':';
    buffer += std::string_view(value);
  } else if constexpr (ConvertibleToBamArrayTag<U>) {
    using V = std::ranges::range_value_t<U>;
    buffer += 'B';
    buffer += ':';
    buffer += ToBamArrayElementType<V>();
    for (const auto& v : value) {
      buffer += ',';
      buffer += std::to_string(v);
    }
  } else if constexpr (std::same_as<U, char>) {
    // A char renders as a single SAM A field holding the character itself. This is handled separately from the other
    // scalars because std::to_string(char) would promote to int and emit the numeric code (e.g. "65" instead of "A").
    buffer += 'A';
    buffer += ':';
    buffer += value;
  } else {
    buffer += ToBamType<U>();
    buffer += ':';
    buffer += std::to_string(value);
  }
}

struct GzipFileDeleter {
  void operator()(gzFile_s* gz) const;
};

using GzipFilePtr = std::unique_ptr<gzFile_s, GzipFileDeleter>;

// Serialize a FASTQ record into `buffer` (cleared first) and write it to `gz`. The caller owns `buffer`
// and should reuse the same instance across calls to avoid per-record allocations; because each worker
// owns its own gzFile and buffer, no synchronization is required (see consensus.cpp).
void GzipWriteFastq(
    gzFile_s* gz, std::string_view read, const vec<u8>& qual, std::string_view read_name, std::string& buffer);

// Serialize a consensus FASTQ record into `buffer` (cleared first) and write it to `gz`. When the consensus aux tags
// are populated, `has_forward`/`has_reverse` select which per-base strand arrays are emitted: the arrays for a strand
// absent from the cluster (all-zero, uninformative) are skipped, while the read-level summary tags are always written.
// A caller that does not know the cluster strandedness should pass true for both to emit every array.
void GzipWriteFastq(gzFile_s* gz,
                    const ConsensusResult& consensus_result,
                    std::string_view read_name,
                    std::string& buffer,
                    bool has_forward,
                    bool has_reverse);

GzipFilePtr OpenGzipFile(const fs::path& path, s32 compression_level = 0);

}  // namespace xoos::read_collapser
