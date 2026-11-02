#include "seq-util.h"

#include <algorithm>

#include "core/variant-type.h"
#include "xoos/types/float.h"
#include "xoos/types/str-container.h"

namespace xoos::svc {

bool IsACTG(const char base) {
  return base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == 'a' || base == 'c' || base == 'g' ||
         base == 't';
}

bool IsNotACTG(const char base) {
  return !IsACTG(base);
}

bool IsAnyNotACTG(const std::string_view seq) {
  return std::ranges::any_of(seq, IsNotACTG);
}

bool ContainsOnlyACTG(const std::string_view seq) {
  return std::ranges::all_of(seq, IsACTG);
}

namespace {
/**
 * @brief Convert a single nucleotide to its 2-bit representation.
 * @param c nucleotide
 * @return 2-bit representation
 */
TwoBit32mer EncodeNucleotide(const char c) {
  switch (c) {
    case 'A':
    case 'a':
      return kTwoBitA;
    case 'C':
    case 'c':
      return kTwoBitC;
    case 'G':
    case 'g':
      return kTwoBitG;
    case 'T':
    case 't':
      return kTwoBitT;
    default:
      return kTwoBitA;
  }
}
}  // namespace

TwoBit32mer SeqToBits(const std::string& seq) {
  static constexpr u64 kMaxSeqLength{32};
  static constexpr u64 kNumBitsPerBase{2};
  TwoBit32mer result{0b00};
  if (seq.length() <= kMaxSeqLength) {
    for (const char n : seq) {
      result <<= kNumBitsPerBase;
      result |= EncodeNucleotide(n);
    }
  } else {
    // longer than 32 nt; use only the first 32 nt
    for (size_t i = 0; i < kMaxSeqLength; ++i) {
      result <<= kNumBitsPerBase;
      result |= EncodeNucleotide(seq[i]);
    }
  }
  return result;
}

f64 SeqToDouble(const std::string& seq) {
  return seq.empty() ? 0 : static_cast<f64>(1 + SeqToBits(seq));
}

bool HasLongerAlt(const std::string_view ref1,
                  const std::string_view alt1,
                  const std::string_view ref2,
                  const std::string_view alt2) {
  // Examples using the REF,ALT representation:
  // `C,CAA` > `C,CA` because `CAA` > `CA`
  // `CAT,C` < `CA,C` because `CA,C` can be formatted as `CAT,CT` and `CT` > `C`
  return alt1.length() > alt2.length() || (alt1.length() == alt2.length() && ref1.length() < ref2.length());
}

std::pair<std::string, std::string> TrimVariant(const std::string& ref, const std::string& alt) {
  const size_t ref_len = ref.length();
  const size_t alt_len = alt.length();
  if (ref_len > 1 && alt_len > 1) {
    size_t count = 0;
    auto ref_itr = ref.rbegin();
    auto alt_itr = alt.rbegin();
    while (ref_itr != ref.rend() && alt_itr != alt.rend()) {
      if (*ref_itr != *alt_itr) {
        break;
      }
      ++ref_itr;
      ++alt_itr;
      ++count;
    }
    return std::make_pair(ref.substr(0, ref_len - count), alt.substr(0, alt_len - count));
  }
  return std::make_pair(ref, alt);
}

namespace {
/**
 * @brief Format two variants that require trimming (SNV+deletion or two indels with shared anchor base).
 * @return Formatted tuple (ref, alt1, alt2) or empty strings if formatting is not possible.
 */
std::tuple<std::string, std::string, std::string> FormatWithTrimming(const std::string& ref1,
                                                                     const std::string& alt1,
                                                                     const std::string& ref2,
                                                                     const std::string& alt2) {
  if (ref1 == ref2) {
    return std::make_tuple(ref1, alt1, alt2);
  }
  // trim both variants to their minimal representation
  const auto& [ref1_trimmed, alt1_trimmed] = TrimVariant(ref1, alt1);
  const auto& [ref2_trimmed, alt2_trimmed] = TrimVariant(ref2, alt2);
  const auto ref1_trimmed_len = ref1_trimmed.length();
  const auto ref2_trimmed_len = ref2_trimmed.length();

  if (ref1_trimmed_len == 1 && ref1_trimmed == ref2_trimmed) {
    // Example
    // insertion 1: A->ATTGG
    // insertion 2: A->ATT
    // reformated:  A->ATTGG,ATT
    return std::make_tuple(ref1_trimmed, alt1_trimmed, alt2_trimmed);
  }
  if (ref1_trimmed_len > ref2_trimmed_len) {
    // Example 1
    // deletion 1:  ATTGG->A
    // deletion 2:  ATT  ->A
    // reformatted: ATTGG->A,AGG

    // Example 2
    // deletion:    ACTGG->A
    // insertion:   A    ->ATTTT
    // reformatted: ACTGG->A,ATTTTCTGG

    // Example 3
    // deletion:    ACTGG->A
    // SNV:         A    ->G
    // reformatted: ACTGG->A,GCTGG
    return std::make_tuple(ref1_trimmed,
                           alt1_trimmed,
                           alt2_trimmed + ref1_trimmed.substr(ref2_trimmed_len, ref1_trimmed_len - ref2_trimmed_len));
  }
  // Example 1
  // deletion 1:  ATT  ->A
  // deletion 2:  ATTGG->A
  // reformatted: ATTGG->AGG,A

  // Example 2
  // insertion:   A    ->ATTTT
  // deletion:    ACTGG->A
  // reformatted: ACTGG->ATTTTCTGG,A

  // Example 3
  // SNV:         A    ->G
  // deletion:    ACTGG->A
  // reformatted: ACTGG->GCTGG,A
  return std::make_tuple(ref2_trimmed,
                         alt1_trimmed + ref2_trimmed.substr(ref1_trimmed_len, ref2_trimmed_len - ref1_trimmed_len),
                         alt2_trimmed);
}

/**
 * @brief Check whether two variant types form a pair of two SNVs, or one SNV and one insertion (either order).
 * @param type1 Variant type of the first variant
 * @param type2 Variant type of the second variant
 * @return true if both types are in {SNV, Insertion} and at least one is an SNV
 */
bool IsSnvOrInsertionPair(const VariantType type1, const VariantType type2) {
  using enum VariantType;
  const bool is_snv_or_ins1 = (type1 == kSNV || type1 == kInsertion);
  const bool is_snv_or_ins2 = (type2 == kSNV || type2 == kInsertion);
  const bool at_least_one_snv = (type1 == kSNV || type2 == kSNV);
  return is_snv_or_ins1 && is_snv_or_ins2 && at_least_one_snv;
}

/**
 * @brief Check whether two variant types form an SNV + deletion pair (either order).
 * @param type1 Variant type of the first variant
 * @param type2 Variant type of the second variant
 * @return true if one type is SNV and the other is Deletion
 */
bool IsSnvDeletionPair(const VariantType type1, const VariantType type2) {
  using enum VariantType;
  return (type1 == kSNV && type2 == kDeletion) || (type1 == kDeletion && type2 == kSNV);
}

/**
 * @brief Check whether two indel variants share the same anchor base (first base of REF and ALT).
 * @param ref1 REF allele for variant 1
 * @param alt1 ALT allele for variant 1
 * @param ref2 REF allele for variant 2
 * @param alt2 ALT allele for variant 2
 * @return true if all four alleles start with the same base
 */
bool ShareAnchorBase(const std::string& ref1,
                     const std::string& alt1,
                     const std::string& ref2,
                     const std::string& alt2) {
  return ref1[0] == alt1[0] && ref2[0] == alt2[0] && ref1[0] == ref2[0];
}
}  // namespace

std::tuple<std::string, std::string, std::string> FormatVariants(const std::string& ref1,
                                                                 const std::string& alt1,
                                                                 const std::string& ref2,
                                                                 const std::string& alt2) {
  if (ref1.empty() || alt1.empty() || ref2.empty() || alt2.empty()) {
    return std::make_tuple("", "", "");
  }
  if (alt1 == "*") {
    return std::make_tuple(ref2, "*", alt2);
  }
  if (alt2 == "*") {
    return std::make_tuple(ref1, alt1, "*");
  }

  const auto type1 = GetVariantType(ref1, alt1);
  const auto type2 = GetVariantType(ref2, alt2);
  if (IsSnvOrInsertionPair(type1, type2) && ref1 == ref2) {
    return std::make_tuple(ref1, alt1, alt2);
  }
  if (IsSnvDeletionPair(type1, type2) || ShareAnchorBase(ref1, alt1, ref2, alt2)) {
    return FormatWithTrimming(ref1, alt1, ref2, alt2);
  }
  // unsupported variant type combination
  return std::make_tuple("", "", "");
}

u32 CountUniqueKmers(const std::string_view seq, const u32 k) {
  StrUnorderedSet unique_kmers;
  const auto num_kmers = static_cast<u32>(seq.length() - k + 1);
  for (u32 i = 0; i < num_kmers; ++i) {
    unique_kmers.emplace(seq.substr(i, k));
  }
  return static_cast<u32>(unique_kmers.size());
}

u32 GetHomopolymerLength(const std::string_view seq) {
  u32 count{0};
  const auto seq_len = static_cast<u32>(seq.length());
  if (seq_len > 0) {
    const char first{seq[0]};
    if (IsNotACTG(first)) {
      return count;
    }
    for (u32 i = 1; i < seq_len; ++i) {
      if (seq[i] == first) {
        ++count;
      } else {
        // not a homopolymer anymore
        break;
      }
    }
    if (count > 0) {
      // increment for the first base
      ++count;
    }
  }
  return count;
}

u32 CountRepeats(const std::string_view seq, const u32 k) {
  u32 count{0};
  const auto seq_len = static_cast<u32>(seq.length());
  if (k < seq_len) {
    const std::string_view repeat = seq.substr(0, k);
    if (IsAnyNotACTG(repeat)) {
      return count;
    }
    u32 i = k;
    while (i + k <= seq_len && repeat == seq.substr(i, k)) {
      ++count;
      i += k;
    }
    if (count > 0) {
      // increment for the first k-mer
      ++count;
    }
  }
  return count;
}

std::optional<u64> FindOverlappingHomopolymer(const std::string& seq,
                                              const u64 pos,
                                              const u64 min_hp_length,
                                              const u64 min_start) {
  constexpr u64 kMinHomopolymerLength = 2;
  // homopolymer must be at least 2 bases long
  if (seq.empty() || min_hp_length < kMinHomopolymerLength || min_hp_length > seq.size()) {
    return std::nullopt;
  }
  if (pos >= seq.size() || pos < min_start) {
    // position out of bounds
    return std::nullopt;
  }

  const char ref{seq[pos]};

  // First look at the bases upstream of the current position
  u64 start = pos;
  if (pos > 0) {
    const u64 lower_bound = std::max(1UL, min_start);
    start = pos - 1;
    while (start >= lower_bound && seq[start] == ref) {
      --start;
    }
    if (start >= lower_bound) {
      // stopped on a mismatch; step forward to the first matching base
      ++start;
    }
  }

  u64 end = pos;
  if (start <= pos) {
    // look at bases downstream of current position
    const u64 upper_bound = std::min(seq.size(), pos + min_hp_length);
    end = pos + 1;
    while (end < upper_bound && seq[end] == ref) {
      ++end;
    }
    if (end < upper_bound) {
      // stopped on a mismatch; step back to the last matching base
      --end;
    }
  }

  if (end - start >= min_hp_length) {
    // homopolymer found
    return start;
  }

  // pos does not overlap a homopolymer
  return std::nullopt;
}

f32 CalculateGcContent(const std::string_view seq) {
  u32 gc_count{0};
  u32 actg_count{0};
  for (const char base : seq) {
    if (IsACTG(base)) {
      ++actg_count;
      if (base == 'G' || base == 'C' || base == 'g' || base == 'c') {
        ++gc_count;
      }
    }
  }
  if (actg_count == 0) {
    return 0.0F;
  }
  return static_cast<f32>(gc_count) / static_cast<f32>(actg_count);
}

}  // namespace xoos::svc
