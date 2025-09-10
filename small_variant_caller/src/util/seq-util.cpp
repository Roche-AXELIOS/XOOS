#include "seq-util.h"

#include <algorithm>

#include "core/variant-type.h"
#include "xoos/types/str-container.h"

namespace xoos::svc {

bool IsACTG(const char base) {
  return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

bool IsNotACTG(const char base) {
  return !IsACTG(base);
}

bool IsAnyNotACTG(const std::string& seq) {
  return std::ranges::any_of(seq, IsNotACTG);
}

bool ContainsOnlyACTG(const std::string& seq) {
  return std::ranges::all_of(seq, IsACTG);
}

/**
 * @brief Convert a single nucleotide to its 2-bit representation.
 * @param c nucleotide
 * @return 2-bit representation
 */
static TwoBit32mer EncodeNucleotide(const char c) {
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

double SeqToDouble(const std::string& seq) {
  return seq.empty() ? 0 : static_cast<double>(1 + SeqToBits(seq));
}

bool HasLongerAlt(const std::string& ref1, const std::string& alt1, const std::string& ref2, const std::string& alt2) {
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

std::tuple<std::string, std::string, std::string> FormatVariants(const std::string& ref1,
                                                                 const std::string& alt1,
                                                                 const std::string& ref2,
                                                                 const std::string& alt2) {
  using enum VariantType;
  if (!ref1.empty() && !alt1.empty() && !ref2.empty() && !alt2.empty()) {
    if (alt1 == "*") {
      return std::make_tuple(ref2, "*", alt2);
    }
    if (alt2 == "*") {
      return std::make_tuple(ref1, alt1, "*");
    }

    const auto type1 = GetVariantType(ref1, alt1);
    const auto type2 = GetVariantType(ref2, alt2);
    if ((type1 == kSNV && type2 == kSNV) || (type1 == kSNV && type2 == kInsertion) ||
        (type1 == kInsertion && type2 == kSNV)) {
      // two SNVs OR SNV + insertion
      if (ref1 == ref2) {
        return std::make_tuple(ref1, alt1, alt2);
      }
    } else if ((type1 == kSNV && type2 == kDeletion) || (type1 == kDeletion && type2 == kSNV) ||
               (ref1[0] == alt1[0] && ref2[0] == alt2[0] && ref1[0] == ref2[0])) {
      // SNV + deletion OR two indels
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
        return std::make_tuple(
            ref1_trimmed,
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
  }
  // cannot be formatted
  return std::make_tuple("", "", "");
}

u32 CountUniqueKmers(const std::string& seq, const u32 k) {
  StrUnorderedSet unique_kmers;
  const auto num_kmers = static_cast<u32>(seq.length() - k + 1);
  for (u32 i = 0; i < num_kmers; ++i) {
    unique_kmers.insert(seq.substr(i, k));
  }
  return static_cast<u32>(unique_kmers.size());
}

u32 GetHomopolymerLength(const std::string& seq) {
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

u32 CountRepeats(const std::string& seq, const u32 k) {
  u32 count{0};
  const auto seq_len = static_cast<u32>(seq.length());
  if (k < seq_len) {
    const std::string& repeat = seq.substr(0, k);
    if (IsAnyNotACTG(repeat)) {
      return count;
    }
    for (u32 i = k; i + k <= seq_len; i += k) {
      if (repeat == seq.substr(i, k)) {
        ++count;
      } else {
        // not a repeat anymore
        break;
      }
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
  // homopolymer must be at least 2 bases long
  if (seq.empty() || min_hp_length < 2 || min_hp_length > seq.size()) {
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
    for (start = pos - 1; start >= std::max(1UL, min_start); --start) {
      if (seq[start] != ref) {
        ++start;
        break;
      }
    }
  }

  u64 end = pos;
  if (start <= pos) {
    // look at bases downstream of current position
    for (end = pos + 1; end < std::min(seq.size(), pos + min_hp_length); ++end) {
      if (seq[end] != ref) {
        --end;
        break;
      }
    }
  }

  if (end - start >= min_hp_length) {
    // homopolymer found
    return start;
  }

  // pos does not overlap a homopolymer
  return std::nullopt;
}

}  // namespace xoos::svc
