#include "phase-engine.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <filesystem>  // NOLINT(build/c++17)
#include <numeric>
#include <optional>
#include <span>
#include <string_view>
#include <unordered_map>
#include <utility>

#include <fmt/ranges.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

#include <taskflow/taskflow.hpp>  // IWYU pragma: keep
#include <taskflow/core/executor.hpp>
#include <taskflow/core/taskflow.hpp>

#include <xoos/io/alignment-reader.h>
#include <xoos/io/htslib-util/htslib-ptr.h>
#include <xoos/io/htslib-util/htslib-util.h>
#include <xoos/log/logging.h>
#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

#include "phaser/params.h"
#include "phaser/phase-block-builder.h"
#include "phaser/variant.h"

namespace xoos::svc::phaser {

// Minimum number of variants required to form a phase group or block (a pair).
constexpr size_t kMinVariantsToPhase = 2;

vec<PhaseGroup> BuildPhaseGroups(const vec<Variant>& variants, const s32 max_distance) {
  vec<PhaseGroup> groups;
  if (variants.empty()) {
    return groups;
  }

  // Sort indices by chrom then position
  vec<s32> sorted_indices(variants.size());
  std::ranges::iota(sorted_indices, 0);
  std::ranges::sort(sorted_indices, [&variants](const s32 a, const s32 b) {
    const auto ua = static_cast<size_t>(a);
    const auto ub = static_cast<size_t>(b);
    if (variants[ua].chrom != variants[ub].chrom) {
      return variants[ua].chrom < variants[ub].chrom;
    }
    return variants[ua].pos < variants[ub].pos;
  });

  // Greedy clustering: extend group while next variant is within max_distance of any variant in group
  size_t i = 0;
  while (i < sorted_indices.size()) {
    PhaseGroup group;
    group.chrom = variants[static_cast<size_t>(sorted_indices[i])].chrom;
    group.variant_indices.push_back(sorted_indices[i]);

    size_t j = i + 1;
    while (j < sorted_indices.size()) {
      const auto& prev = variants[static_cast<size_t>(group.variant_indices.back())];
      const auto& curr = variants[static_cast<size_t>(sorted_indices[j])];
      if (curr.chrom != prev.chrom || curr.pos - prev.pos > max_distance) {
        break;
      }
      group.variant_indices.push_back(sorted_indices[j]);
      ++j;
    }

    // Only keep groups with 2+ variants
    if (group.variant_indices.size() >= kMinVariantsToPhase) {
      groups.push_back(std::move(group));
    }
    i = j;
  }

  Logging::Info("Built {} candidate phase groups from {} variants", groups.size(), variants.size());
  return groups;
}

namespace {

// Aligned pair: (query_offset, ref_pos). -1 means gap.
struct AlignedPair {
  // -1 if deletion
  s32 query_pos;
  // -1 if insertion
  s64 ref_pos;
};

// A contiguous run of anchor positions in the aligned-pair list, starting at `start` for `len` bases.
struct AnchorRun {
  s32 start;
  s32 len;
};

// Emit aligned pairs for a CIGAR match run (both reference and query advance).
void AppendMatchRun(vec<AlignedPair>& pairs, const s32 len, s32& seq_pos, s64& ref_pos) {
  for (s32 k = 0; k < len; ++k) {
    pairs.push_back({seq_pos, ref_pos});
    ++seq_pos;
    ++ref_pos;
  }
}

// Emit aligned pairs for an insertion run (only the query advances).
void AppendInsertionRun(vec<AlignedPair>& pairs, const s32 len, s32& seq_pos) {
  for (s32 k = 0; k < len; ++k) {
    pairs.push_back({seq_pos, -1});
    ++seq_pos;
  }
}

// Emit aligned pairs for a deletion/skip run (only the reference advances).
void AppendDeletionRun(vec<AlignedPair>& pairs, const s32 len, s64& ref_pos) {
  for (s32 k = 0; k < len; ++k) {
    pairs.push_back({-1, ref_pos});
    ++ref_pos;
  }
}

vec<AlignedPair> BuildAlignedPairs(const std::span<const u32> cigar, const s64 read_pos) {
  vec<AlignedPair> pairs;
  s64 ref_pos = read_pos;
  s32 seq_pos = 0;

  for (const u32 cigar_op : cigar) {
    const auto op = static_cast<s32>(bam_cigar_op(cigar_op));
    const auto len = static_cast<s32>(bam_cigar_oplen(cigar_op));

    switch (op) {
      case BAM_CMATCH:
      case BAM_CEQUAL:
      case BAM_CDIFF:
        AppendMatchRun(pairs, len, seq_pos, ref_pos);
        break;
      case BAM_CINS:
        AppendInsertionRun(pairs, len, seq_pos);
        break;
      case BAM_CDEL:
      case BAM_CREF_SKIP:
        AppendDeletionRun(pairs, len, ref_pos);
        break;
      case BAM_CSOFT_CLIP:
        seq_pos += len;
        break;
      case BAM_CHARD_CLIP:
      case BAM_CPAD:
      default:
        break;
    }
  }
  return pairs;
}

// Check base quality at all bases involved in a variant call.
// For SNVs: checks the single substituted base.
// For MNVs: checks all substituted bases.
// For insertions: checks anchor base(s) AND inserted bases.
// For deletions: checks anchor base(s) flanking the deletion.
// Returns true if ANY variant base is below min_bq (i.e., BQ fails).
// Locate the index in `pairs` whose reference position equals `var_pos`.
// Returns -1 when no aligned pair covers the position.
s32 FindAnchorStart(const vec<AlignedPair>& pairs, const s64 var_pos) {
  for (s32 i = 0; i < static_cast<s32>(pairs.size()); ++i) {
    if (pairs[static_cast<size_t>(i)].ref_pos == var_pos) {
      return i;
    }
  }
  return -1;
}

// SNV base quality check: the single base at var_pos must meet min_bq.
// Returns true when the read fails the check (low BQ or position absent).
bool SnvBQFails(const vec<AlignedPair>& pairs, const std::span<const u8> qual, const s64 var_pos, const s32 min_bq) {
  for (const auto& p : pairs) {
    if (p.ref_pos == var_pos && p.query_pos >= 0) {
      return std::cmp_less(qual[static_cast<size_t>(p.query_pos)], min_bq);
    }
  }
  // Position not found → fail.
  return true;
}

// MNV base quality check: every base across the span must meet min_bq.
bool MnvBQFails(const vec<AlignedPair>& pairs,
                const std::span<const u8> qual,
                const s64 var_pos,
                const size_t span,
                const s32 min_bq) {
  for (size_t b = 0; b < span; ++b) {
    const s64 target = var_pos + static_cast<s64>(b);
    bool found = false;
    for (const auto& p : pairs) {
      if (p.ref_pos == target && p.query_pos >= 0) {
        if (std::cmp_less(qual[static_cast<size_t>(p.query_pos)], min_bq)) {
          return true;
        }
        found = true;
        break;
      }
    }
    if (!found) {
      return true;
    }
  }
  return false;
}

// Check a run of anchor bases (`anchor.start` for `anchor.len` positions).
// Returns true when any base fails (out of range, unaligned, or low BQ).
bool AnchorBasesFail(const vec<AlignedPair>& pairs,
                     const std::span<const u8> qual,
                     const AnchorRun anchor,
                     const s32 min_bq) {
  for (s32 k = 0; k < anchor.len; ++k) {
    const s32 idx = anchor.start + k;
    if (idx >= static_cast<s32>(pairs.size()) || pairs[static_cast<size_t>(idx)].query_pos < 0) {
      return true;
    }
    if (std::cmp_less(qual[static_cast<size_t>(pairs[static_cast<size_t>(idx)].query_pos)], min_bq)) {
      return true;
    }
  }
  return false;
}

// Insertion base quality check: anchor bases plus the inserted bases must meet min_bq.
bool InsertionBQFails(const vec<AlignedPair>& pairs,
                      const std::span<const u8> qual,
                      const s64 var_pos,
                      const size_t ref_len,
                      const size_t alt_len,
                      const s32 min_bq) {
  const auto anchor_len = static_cast<s32>(ref_len);
  const auto ins_len = static_cast<s32>(alt_len - ref_len);

  const s32 anchor_start = FindAnchorStart(pairs, var_pos);
  if (anchor_start < 0 || AnchorBasesFail(pairs, qual, {.start = anchor_start, .len = anchor_len}, min_bq)) {
    return true;
  }

  for (s32 k = 0; k < ins_len; ++k) {
    const s32 idx = anchor_start + anchor_len + k;
    if (idx >= static_cast<s32>(pairs.size())) {
      return true;
    }
    const auto& pair = pairs[static_cast<size_t>(idx)];
    if (pair.query_pos < 0 || pair.ref_pos >= 0) {
      // No insertion here.
      break;
    }
    if (std::cmp_less(qual[static_cast<size_t>(pair.query_pos)], min_bq)) {
      return true;
    }
  }
  return false;
}

// Deletion base quality check: the anchor bases must meet min_bq.
bool DeletionBQFails(const vec<AlignedPair>& pairs,
                     const std::span<const u8> qual,
                     const s64 var_pos,
                     const size_t alt_len,
                     const s32 min_bq) {
  const s32 anchor_start = FindAnchorStart(pairs, var_pos);
  if (anchor_start < 0) {
    return true;
  }
  return AnchorBasesFail(pairs, qual, {.start = anchor_start, .len = static_cast<s32>(alt_len)}, min_bq);
}

bool CheckVariantBQ(const std::span<const u8> qual,
                    const std::span<const u32> cigar,
                    const s64 read_pos,
                    const Variant& var,
                    const s32 min_bq) {
  const vec<AlignedPair> pairs = BuildAlignedPairs(cigar, read_pos);

  if (var.IsSnv()) {
    return SnvBQFails(pairs, qual, var.pos, min_bq);
  }
  if (var.IsMnv()) {
    return MnvBQFails(pairs, qual, var.pos, var.ref.size(), min_bq);
  }
  if (var.IsInsertion()) {
    return InsertionBQFails(pairs, qual, var.pos, var.ref.size(), var.alt.size(), min_bq);
  }
  if (var.IsDeletion()) {
    return DeletionBQFails(pairs, qual, var.pos, var.alt.size(), min_bq);
  }
  // unknown variant type → fail conservatively
  return true;
}

// Return the index into `pairs` whose reference position equals `var_pos`, or -1 if none.
s32 FindPairAtRefPos(const vec<AlignedPair>& pairs, const s64 var_pos) {
  for (s32 i = 0; i < static_cast<s32>(pairs.size()); ++i) {
    if (pairs[static_cast<size_t>(i)].ref_pos == var_pos) {
      return i;
    }
  }
  return -1;
}

// SNV allele call.
Allele DetermineSnvAllele(const u8* const seq, const vec<AlignedPair>& pairs, const Variant& var) {
  using enum Allele;
  const char ref = var.ref[0];
  const char alt = var.alt[0];
  for (const auto& p : pairs) {
    if (p.ref_pos == var.pos) {
      if (p.query_pos < 0) {
        return kUnknown;
      }
      const char base = seq_nt16_str[bam_seqi(seq, p.query_pos)];
      if (base == alt) {
        return kAlt;
      }
      if (base == ref) {
        return kRef;
      }
      return kUnknown;
    }
  }
  return kUnknown;
}

// MNV allele call (ref and alt equal length, >1 base).
Allele DetermineMnvAllele(const u8* const seq, const vec<AlignedPair>& pairs, const Variant& var) {
  using enum Allele;
  const std::string_view ref = var.ref;
  const std::string_view alt = var.alt;
  std::unordered_map<s64, s32> ref_to_query;
  for (const auto& p : pairs) {
    if (p.ref_pos >= var.pos && p.ref_pos < var.pos + static_cast<s64>(ref.size())) {
      ref_to_query[p.ref_pos] = p.query_pos;
    }
  }
  if (ref_to_query.size() < ref.size()) {
    return kUnknown;
  }
  bool matches_alt = true;
  bool matches_ref = true;
  for (size_t b = 0; b < ref.size(); ++b) {
    const auto it = ref_to_query.find(var.pos + static_cast<s64>(b));
    if (it == ref_to_query.end() || it->second < 0) {
      return kUnknown;
    }
    const char base = seq_nt16_str[bam_seqi(seq, it->second)];
    if (base != alt[b]) {
      matches_alt = false;
    }
    if (base != ref[b]) {
      matches_ref = false;
    }
  }
  if (matches_alt) {
    return kAlt;
  }
  if (matches_ref) {
    return kRef;
  }
  return kUnknown;
}

// Collect `count` anchor bases starting at `anchor_start`; returns nullopt if any base is unaligned/out of range.
std::optional<std::string> CollectAnchorBases(const u8* const seq,
                                              const vec<AlignedPair>& pairs,
                                              const s32 anchor_start,
                                              const s32 count,
                                              const bool require_aligned_ref) {
  std::string bases;
  for (s32 k = 0; k < count; ++k) {
    const s32 idx = anchor_start + k;
    if (idx >= static_cast<s32>(pairs.size()) || pairs[static_cast<size_t>(idx)].query_pos < 0) {
      return std::nullopt;
    }
    if (require_aligned_ref && pairs[static_cast<size_t>(idx)].ref_pos < 0) {
      return std::nullopt;
    }
    bases += seq_nt16_str[bam_seqi(seq, pairs[static_cast<size_t>(idx)].query_pos)];
  }
  return bases;
}

// Deletion allele call (ref longer than alt, e.g. ref=TA alt=T).
Allele DetermineDeletionAllele(const u8* const seq, const vec<AlignedPair>& pairs, const Variant& var) {
  using enum Allele;
  const std::string_view ref = var.ref;
  const std::string_view alt = var.alt;
  const auto anchor_len = static_cast<s32>(alt.size());
  const auto del_len = static_cast<s32>(ref.size() - alt.size());

  const s32 anchor_start = FindPairAtRefPos(pairs, var.pos);
  if (anchor_start < 0) {
    return kUnknown;
  }

  const std::optional<std::string> anchor_bases = CollectAnchorBases(seq, pairs, anchor_start, anchor_len, true);
  if (!anchor_bases) {
    return kUnknown;
  }

  // Check for deletion: next del_len pairs should be (query=-1, ref=pos).
  s32 del_count = 0;
  for (s32 k = 0; k < del_len; ++k) {
    const s32 idx = anchor_start + anchor_len + k;
    const bool is_deletion_pair = idx < static_cast<s32>(pairs.size()) &&
                                  pairs[static_cast<size_t>(idx)].query_pos < 0 &&
                                  pairs[static_cast<size_t>(idx)].ref_pos >= 0;
    if (!is_deletion_pair) {
      break;
    }
    ++del_count;
  }

  if (*anchor_bases == alt && del_count == del_len) {
    return kAlt;
  }

  // Check for REF: all ref bases present.
  const std::optional<std::string> ref_bases =
      CollectAnchorBases(seq, pairs, anchor_start, static_cast<s32>(ref.size()), false);
  if (ref_bases && *ref_bases == ref) {
    return kRef;
  }

  return kUnknown;
}

// Insertion allele call (alt longer than ref, e.g. ref=C alt=CT).
Allele DetermineInsertionAllele(const u8* const seq, const vec<AlignedPair>& pairs, const Variant& var) {
  using enum Allele;
  const std::string_view ref = var.ref;
  const std::string_view alt = var.alt;
  const auto anchor_len = static_cast<s32>(ref.size());
  const auto ins_len = static_cast<s32>(alt.size() - ref.size());

  const s32 anchor_start = FindPairAtRefPos(pairs, var.pos);
  if (anchor_start < 0) {
    return kUnknown;
  }

  const std::optional<std::string> anchor_bases = CollectAnchorBases(seq, pairs, anchor_start, anchor_len, false);
  if (!anchor_bases || *anchor_bases != ref.substr(0, static_cast<size_t>(anchor_len))) {
    return kUnknown;
  }

  // Check for insertion: next ins_len pairs should be (query=pos, ref=-1).
  std::string ins_bases;
  for (s32 k = 0; k < ins_len; ++k) {
    const s32 idx = anchor_start + anchor_len + k;
    const bool is_insertion_pair = idx < static_cast<s32>(pairs.size()) &&
                                   pairs[static_cast<size_t>(idx)].query_pos >= 0 &&
                                   pairs[static_cast<size_t>(idx)].ref_pos < 0;
    if (!is_insertion_pair) {
      break;
    }
    ins_bases += seq_nt16_str[bam_seqi(seq, pairs[static_cast<size_t>(idx)].query_pos)];
  }

  if (static_cast<s32>(ins_bases.size()) == ins_len && (*anchor_bases + ins_bases) == alt) {
    return kAlt;
  }

  // No insertion → REF.
  const s32 next_idx = anchor_start + anchor_len;
  if (next_idx < static_cast<s32>(pairs.size()) && pairs[static_cast<size_t>(next_idx)].ref_pos >= 0) {
    return kRef;
  }

  return kUnknown;
}

}  // namespace

Allele DetermineAllele(
    const u8* const seq, const u32* const cigar, const s32 n_cigar, const s64 read_pos, const Variant& var) {
  using enum Allele;
  const vec<AlignedPair> pairs = BuildAlignedPairs({cigar, static_cast<size_t>(n_cigar)}, read_pos);

  if (var.IsSnv()) {
    return DetermineSnvAllele(seq, pairs, var);
  }
  if (var.IsMnv()) {
    return DetermineMnvAllele(seq, pairs, var);
  }
  if (var.IsDeletion()) {
    return DetermineDeletionAllele(seq, pairs, var);
  }
  if (var.IsInsertion()) {
    return DetermineInsertionAllele(seq, pairs, var);
  }

  return kUnknown;
}

bool VariantsOverlap(const Variant& a, const Variant& b) {
  if (a.chrom != b.chrom) {
    return false;
  }
  // Check if either variant's ref span covers the other's position.
  // NOSONAR(cpp:M23_002) — symmetric half-open interval overlap; both clauses are
  // required for correctness as a public predicate. The first clause reads as
  // always-true only because the internal caller passes position-sorted variants.
  return (a.pos < b.EndPos() && b.pos < a.EndPos());
}

namespace {

// Extract the read's sequence across a reference span, including insertions and accounting for deletions.
std::string ExtractReadSequence(const vec<AlignedPair>& pairs,
                                const u8* const seq,
                                const s64 span_start,
                                const s64 span_end) {
  std::string result;
  bool in_span = false;

  for (const auto& p : pairs) {
    if (p.ref_pos >= 0) {
      if (p.ref_pos >= span_start && p.ref_pos < span_end) {
        in_span = true;
        if (p.query_pos >= 0) {
          result += seq_nt16_str[bam_seqi(seq, p.query_pos)];
        }
        // Deletion (query_pos < 0): base absent from read, skip.
      } else if (p.ref_pos >= span_end) {
        break;
      }
    } else if (in_span && p.query_pos >= 0) {
      // Insertion within the span (ref_pos < 0 in this branch).
      result += seq_nt16_str[bam_seqi(seq, p.query_pos)];
    }
  }

  return result;
}

// Build expected haplotype for overlapping variant pair.
// Constructs the sequence that a read should show for a given combination of REF/ALT alleles.
//
// For overlapping variants where first's REF span covers second's position:
// - ALT/ALT: first's deletion removes second's anchor base, so second's ALT contribution
//   is only the inserted bases (ALT minus anchor). The aligner represents this as
//   deletion + insertion + match at the next reference position.
// - ALT/REF: first's deletion removes second's position entirely. Physically contradictory
//   but included for completeness (won't match real reads).
// - REF/ALT: first stays REF, second's ALT is applied normally at its position.
//
// Parameters:
//   ref_seq: reference sequence covering [span_start, span_start + ref_seq.size())
//            includes 1 context base beyond the variants' reference spans
//   span_start: start of the reference span
//   first/second: variants ordered by position (first->pos <= second->pos)
//   first_allele/second_allele: which allele (REF or ALT) to apply for each variant
// Build the haplotype when `second` starts within `first`'s REF span (true overlap).
std::string BuildOverlappingHaplotype(const std::string_view ref_seq,
                                      const s64 span_start,
                                      const Variant* const first,
                                      const Variant* const second,
                                      const Allele first_allele,
                                      const Allele second_allele) {
  using enum Allele;
  const bool apply_first = first_allele == kAlt;
  const bool apply_second = second_allele == kAlt;
  // Both REF: the reference sequence is unchanged.
  if (!apply_first && !apply_second) {
    return std::string{ref_seq};
  }

  std::string result;
  if (apply_first && apply_second) {
    // Both ALT: first's deletion removes second's anchor base.
    // Combined haplotype: ref_before_first + first_ALT + second_inserted_bases + ref_after,
    // where second_inserted_bases = second->alt minus the anchor (second->ref).
    const s64 first_offset = first->pos - span_start;
    result = ref_seq.substr(0, static_cast<size_t>(first_offset));
    result += first->alt;
    // For an insertion like G→GTTTT the anchor G is deleted by first, so only TTTT is inserted.
    if (second->alt.size() > second->ref.size()) {
      result += second->alt.substr(second->ref.size());
    }
    // ref_seq spans both variants plus one trailing context base, so `after`
    // never exceeds ref_seq.size(); substr safely returns empty at the boundary.
    const s64 after = std::max(first->EndPos(), second->EndPos()) - span_start;
    result += ref_seq.substr(static_cast<size_t>(after));
  } else if (apply_first) {
    // First ALT only: replaces first's entire REF span (which includes second's position).
    // Physically contradictory for true overlaps — won't match real reads.
    const s64 first_offset = first->pos - span_start;
    result = ref_seq.substr(0, static_cast<size_t>(first_offset));
    result += first->alt;
    const s64 after_first = first->EndPos() - span_start;
    if (after_first < static_cast<s64>(ref_seq.size())) {
      result += ref_seq.substr(static_cast<size_t>(after_first));
    }
  } else {
    // Second ALT only: first stays REF, apply second normally.
    const s64 second_offset = second->pos - span_start;
    result = ref_seq.substr(0, static_cast<size_t>(second_offset));
    result += second->alt;
    // ref_seq spans both variants plus one trailing context base, so `after`
    // never exceeds ref_seq.size(); substr safely returns empty at the boundary.
    const s64 after = std::max(first->EndPos(), second->EndPos()) - span_start;
    result += ref_seq.substr(static_cast<size_t>(after));
  }
  return result;
}

// Build the haplotype when the two variants do not truly overlap — sequential substitution.
std::string BuildSequentialHaplotype(const std::string_view ref_seq,
                                     const s64 span_start,
                                     const Variant* const first,
                                     const Variant* const second,
                                     const Allele first_allele,
                                     const Allele second_allele) {
  using enum Allele;
  const bool apply_first = first_allele == kAlt;
  const bool apply_second = second_allele == kAlt;
  std::string result;
  s64 cursor = span_start;
  if (first->pos > cursor) {
    result = ref_seq.substr(0, static_cast<size_t>(first->pos - span_start));
  }
  result += apply_first ? first->alt : first->ref;
  cursor = first->EndPos();
  if (second->pos > cursor) {
    result += ref_seq.substr(static_cast<size_t>(cursor - span_start), static_cast<size_t>(second->pos - cursor));
  }
  result += apply_second ? second->alt : second->ref;
  cursor = second->EndPos();
  if (cursor - span_start < static_cast<s64>(ref_seq.size())) {
    result += ref_seq.substr(static_cast<size_t>(cursor - span_start));
  }
  return result;
}

std::string BuildOverlapHaplotype(const std::string_view ref_seq,
                                  const s64 span_start,
                                  const Variant* const first,
                                  const Variant* const second,
                                  const Allele first_allele,
                                  const Allele second_allele) {
  if (second->pos < first->EndPos()) {
    return BuildOverlappingHaplotype(ref_seq, span_start, first, second, first_allele, second_allele);
  }
  return BuildSequentialHaplotype(ref_seq, span_start, first, second, first_allele, second_allele);
}

// Look up the read base aligned to reference position `ref_pos`.
// Returns std::nullopt when no aligned query base covers the position.
std::optional<char> ReadBaseAt(const vec<AlignedPair>& pairs, const u8* const seq, const s64 ref_pos) {
  for (const auto& p : pairs) {
    if (p.ref_pos == ref_pos && p.query_pos >= 0) {
      return seq_nt16_str[bam_seqi(seq, p.query_pos)];
    }
  }
  return std::nullopt;
}

// Reconstruct the reference sequence across the span covering both ordered variants,
// filling any gap between them and appending one trailing context base from the read.
// Returns std::nullopt when a required read base is missing.
std::optional<std::string> BuildOverlapRefSeq(const vec<AlignedPair>& pairs,
                                              const u8* const seq,
                                              const Variant* const first,
                                              const Variant* const second,
                                              const s64 span_end) {
  std::string ref_seq = first->ref;
  if (second->EndPos() > first->EndPos()) {
    if (second->pos >= first->EndPos()) {
      for (s64 rp = first->EndPos(); rp < second->pos; ++rp) {
        const auto base = ReadBaseAt(pairs, seq, rp);
        if (!base.has_value()) {
          return std::nullopt;
        }
        ref_seq += *base;
      }
      ref_seq += second->ref;
    } else {
      const s64 extra_start = first->EndPos() - second->pos;
      if (extra_start < static_cast<s64>(second->ref.size())) {
        ref_seq += second->ref.substr(static_cast<size_t>(extra_start));
      }
    }
  }

  // Append context base from the read at span_end - 1.
  const auto context_base = ReadBaseAt(pairs, seq, span_end - 1);
  if (!context_base.has_value()) {
    return std::nullopt;
  }
  ref_seq += *context_base;
  return ref_seq;
}

// Match the read sequence against the four haplotypes and return the (a, b) allele pair.
std::pair<Allele, Allele> MatchOverlapHaplotype(const std::string& read_seq,
                                                const std::string& ref_seq,
                                                const s64 span_start,
                                                const Variant* const first,
                                                const Variant* const second,
                                                const bool first_is_a) {
  using enum Allele;
  const std::string hap_ref_ref = BuildOverlapHaplotype(ref_seq, span_start, first, second, kRef, kRef);
  const std::string hap_alt_alt = BuildOverlapHaplotype(ref_seq, span_start, first, second, kAlt, kAlt);
  const std::string hap_alt_ref = BuildOverlapHaplotype(ref_seq, span_start, first, second, kAlt, kRef);
  const std::string hap_ref_alt = BuildOverlapHaplotype(ref_seq, span_start, first, second, kRef, kAlt);

  Logging::Debug("Overlap haplotypes: ref_ref='{}' alt_alt='{}' alt_ref='{}' ref_alt='{}' read='{}'",
                 hap_ref_ref,
                 hap_alt_alt,
                 hap_alt_ref,
                 hap_ref_alt,
                 read_seq);

  if (read_seq == hap_alt_alt) {
    return {kAlt, kAlt};
  }
  if (read_seq == hap_ref_ref) {
    return {kRef, kRef};
  }
  if (read_seq == hap_alt_ref) {
    return first_is_a ? std::make_pair(kAlt, kRef) : std::make_pair(kRef, kAlt);
  }
  if (read_seq == hap_ref_alt) {
    return first_is_a ? std::make_pair(kRef, kAlt) : std::make_pair(kAlt, kRef);
  }
  return {kUnknown, kUnknown};
}

}  // namespace

std::pair<Allele, Allele> DetermineAllelesOverlapping(const u8* const seq,
                                                      const u32* const cigar,
                                                      const s32 n_cigar,
                                                      const s64 read_pos,
                                                      const Variant& var_a,
                                                      const Variant& var_b) {
  using enum Allele;
  const vec<AlignedPair> pairs = BuildAlignedPairs({cigar, static_cast<size_t>(n_cigar)}, read_pos);

  // Order variants by position.
  const Variant* first = &var_a;
  const Variant* second = &var_b;
  if (var_b.pos < var_a.pos) {
    std::swap(first, second);
  }

  // Define the reference span covering both variants, plus 1 context base.
  // The context base is needed because aligners may represent complex indel events
  // (e.g., 3D4I instead of 3D5I) where the last alt base aligns as a match at the
  // next reference position rather than as part of the insertion.
  const s64 span_start = first->pos;
  const s64 span_end = std::max(first->EndPos(), second->EndPos()) + 1;

  const auto ref_seq = BuildOverlapRefSeq(pairs, seq, first, second, span_end);
  if (!ref_seq.has_value()) {
    return {kUnknown, kUnknown};
  }

  const std::string read_seq = ExtractReadSequence(pairs, seq, span_start, span_end);
  return MatchOverlapHaplotype(read_seq, *ref_seq, span_start, first, second, first == &var_a);
}

namespace {

// Per-read accumulator for split-read merging.
// Mutect2 bamout hard-clips reads at active region boundaries, so a single
// physical read may appear as multiple BAM records. We group by read name
// and keep the best allele call from any record covering each variant.
//
// Mutect2 bamout also emits the same read aligned to different haplotypes,
// so the same read name can appear with different allele calls. We prefer
// calls from records that cover both variants (joint haplotype context)
// over calls from records covering only one variant.
struct ReadAlleleInfo {
  // best allele call for var_a (kUnknown = undetermined)
  Allele allele_a = Allele::kUnknown;
  // best allele call for var_b (kUnknown = undetermined)
  Allele allele_b = Allele::kUnknown;
  // BQ too low at var_a
  bool bq_fail_a = false;
  // BQ too low at var_b
  bool bq_fail_b = false;
  // allele_a came from a record covering both variants
  bool from_joint_a = false;
  // allele_b came from a record covering both variants
  bool from_joint_b = false;
};

// Allele calls plus base-quality outcome for one BAM record against a variant pair.
struct RecordAlleleCall {
  Allele allele_a = Allele::kUnknown;
  Allele allele_b = Allele::kUnknown;
  bool bq_fail_a = false;
  bool bq_fail_b = false;
  bool covers_a = false;
  bool covers_b = false;
};

// Return true if the record should be skipped (not a usable tumor read).
bool ShouldSkipRecord(const bam1_t* const aln,
                      const std::string_view tumor_rg,
                      const std::string_view tumor_rg_prefix) {
  // htslib BAM_F* flag constants are signed int macros; cast to unsigned so the mask is a
  // clean unsigned bitwise op (matches compute-bam-region-features.cpp).
  static constexpr auto kSkipMask = static_cast<u32>(BAM_FUNMAP) | static_cast<u32>(BAM_FSECONDARY) |
                                    static_cast<u32>(BAM_FSUPPLEMENTARY) | static_cast<u32>(BAM_FDUP) |
                                    static_cast<u32>(BAM_FQCFAIL);
  if ((static_cast<u32>(aln->core.flag) & kSkipMask) != 0) {
    return true;
  }
  // Filter by read group — only tumor reads. Use prefix match to handle
  // Mutect2 bamout RG suffixes (e.g. HPSID_007-0C134A2E).
  const u8* const rg_tag = bam_aux_get(aln, "RG");
  if (rg_tag == nullptr) {
    return true;
  }
  // A single-sample VCF can omit --tumor-sample-name; in that mode all tagged
  // BAM reads belong to the tumor sample and are eligible for phasing.
  if (tumor_rg.empty()) {
    return false;
  }
  const std::string_view rg_view(bam_aux2Z(rg_tag));
  return rg_view != tumor_rg && !rg_view.starts_with(tumor_rg_prefix);
}

// Determine the allele calls (and base-quality outcome) for one BAM record.
RecordAlleleCall CallRecordAlleles(bam1_t* const aln,
                                   const Variant& var_a,
                                   const Variant& var_b,
                                   const bool use_overlap_mode,
                                   const s32 min_base_quality) {
  using enum Allele;
  RecordAlleleCall call;
  const s64 read_start = aln->core.pos;
  const s64 read_end = bam_endpos(aln);

  call.covers_a = (read_start <= var_a.pos && read_end > var_a.pos);
  call.covers_b = (read_start <= var_b.pos && read_end > var_b.pos);
  if (!call.covers_a && !call.covers_b) {
    return call;
  }

  const u8* const seq = bam_get_seq(aln);
  const u32* const cigar = bam_get_cigar(aln);
  const auto n_cigar = static_cast<s32>(aln->core.n_cigar);

  if (call.covers_a && call.covers_b && use_overlap_mode) {
    const auto [a, b] = DetermineAllelesOverlapping(seq, cigar, n_cigar, read_start, var_a, var_b);
    call.allele_a = a;
    call.allele_b = b;
  } else {
    if (call.covers_a) {
      call.allele_a = DetermineAllele(seq, cigar, n_cigar, read_start, var_a);
    }
    if (call.covers_b) {
      call.allele_b = DetermineAllele(seq, cigar, n_cigar, read_start, var_b);
    }
  }

  // Check base quality at all bases involved in each variant.
  // For indels this includes inserted bases, not just the anchor position.
  if (min_base_quality > 0) {
    const std::span<const u8> qual{bam_get_qual(aln), static_cast<size_t>(aln->core.l_qseq)};
    const std::span<const u32> cigar_span{cigar, static_cast<size_t>(n_cigar)};
    // allele_a/allele_b are only determined when the record covers the respective
    // variant, so a determined allele already implies coverage.
    if (call.allele_a != kUnknown) {
      call.bq_fail_a = CheckVariantBQ(qual, cigar_span, read_start, var_a, min_base_quality);
    }
    if (call.allele_b != kUnknown) {
      call.bq_fail_b = CheckVariantBQ(qual, cigar_span, read_start, var_b, min_base_quality);
    }
  }
  return call;
}

// Merge one record's allele calls into the per-read accumulator.
// Prefer calls from records covering both variants (joint context) over
// single-variant records, since Mutect2 bamout can emit the same read
// aligned to different haplotypes.
void MergeRecordCall(ReadAlleleInfo& info, const RecordAlleleCall& call) {
  using enum Allele;
  const bool is_joint = call.covers_a && call.covers_b;

  if (call.allele_a != kUnknown && !call.bq_fail_a) {
    if (info.allele_a == kUnknown || (is_joint && !info.from_joint_a)) {
      info.allele_a = call.allele_a;
      info.bq_fail_a = false;
      info.from_joint_a = is_joint;
    }
  } else if (call.allele_a != kUnknown && info.allele_a == kUnknown) {
    // Reaching here implies call.bq_fail_a is set (the branch above requires !bq_fail).
    info.bq_fail_a = true;
  }

  if (call.allele_b != kUnknown && !call.bq_fail_b) {
    if (info.allele_b == kUnknown || (is_joint && !info.from_joint_b)) {
      info.allele_b = call.allele_b;
      info.bq_fail_b = false;
      info.from_joint_b = is_joint;
    }
  } else if (call.allele_b != kUnknown && info.allele_b == kUnknown) {
    // Reaching here implies call.bq_fail_b is set (the branch above requires !bq_fail).
    info.bq_fail_b = true;
  }
}

// Tally cis/trans/ref counts from accumulated per-read allele calls into result.
void TallyPhasingCounts(const StrUnorderedMap<ReadAlleleInfo>& read_alleles, PhasingResult& result) {
  using enum Allele;
  for (const auto& [qname, info] : read_alleles) {
    if (info.allele_a == kUnknown || info.allele_b == kUnknown || info.bq_fail_a || info.bq_fail_b) {
      continue;
    }
    ++result.total_spanning;
    if (info.allele_a == kAlt && info.allele_b == kAlt) {
      ++result.cis_count;
    } else if (info.allele_a == kRef && info.allele_b == kRef) {
      ++result.ref_count;
    } else {
      ++result.trans_count;
    }
  }
}

}  // namespace

PhasingResult PhaseVariantPair(const io::AlignmentReader& bam_handle,
                               const std::string& tumor_rg,
                               const Variant& var_a,
                               const Variant& var_b,
                               const PhasingThresholds& thresholds) {
  PhasingResult result;

  if (bam_handle.bam_file == nullptr || bam_handle.hdr == nullptr || bam_handle.idx == nullptr) {
    Logging::Warn("BAM handle not open");
    return result;
  }

  // Pre-compute the RG prefix for matching (avoid repeated string concat in hot loop)
  const std::string tumor_rg_prefix = tumor_rg + "-";

  // Region spanning both variants
  const s64 region_start = std::min(var_a.pos, var_b.pos);
  const s64 region_end = std::max(var_a.EndPos(), var_b.EndPos());
  const s32 tid = sam_hdr_name2tid(bam_handle.hdr.get(), var_a.chrom.c_str());
  if (tid < 0) {
    Logging::Warn("Contig {} not found in BAM", var_a.chrom);
    return result;
  }

  // unique_ptr with custom deleter adopts the htslib C-API pointer; make_unique cannot be used.
  const io::HtsItrPtr iter{// NOSONAR(cpp:M23_286)
                           sam_itr_queryi(bam_handle.idx.get(), tid, region_start, region_end)};
  if (iter == nullptr) {
    return result;
  }

  // Use haplotype-based comparison when variants overlap in reference coordinates
  const bool use_overlap_mode = VariantsOverlap(var_a, var_b);
  if (use_overlap_mode) {
    Logging::Debug("Using overlap mode for {}:{} and {}:{}", var_a.chrom, var_a.pos + 1, var_b.chrom, var_b.pos + 1);
  }

  // Accumulate allele calls per read name to handle split reads.
  // Each BAM record contributes allele info for whichever variant(s) it covers.
  StrUnorderedMap<ReadAlleleInfo> read_alleles;

  const io::Bam1Ptr aln_ptr = io::BamInit1();
  bam1_t* const aln = aln_ptr.get();
  while (sam_itr_next(bam_handle.bam_file.get(), iter.get(), aln) >= 0) {
    if (ShouldSkipRecord(aln, tumor_rg, tumor_rg_prefix)) {
      continue;
    }

    const RecordAlleleCall call = CallRecordAlleles(aln, var_a, var_b, use_overlap_mode, thresholds.min_base_quality);
    if (!call.covers_a && !call.covers_b) {
      continue;
    }

    const char* const qname = bam_get_qname(aln);
    MergeRecordCall(read_alleles[qname], call);
  }

  // Count spanning reads from accumulated per-read allele calls
  TallyPhasingCounts(read_alleles, result);

  // Apply phasing criteria — phase as cis or trans
  const s32 informative = result.cis_count + result.trans_count;
  if (result.total_spanning >= thresholds.min_shared_reads && informative > 0) {
    if (result.cis_count >= thresholds.min_cis_reads && result.CisFraction() >= thresholds.min_cis_fraction) {
      result.is_phased = true;
      result.is_cis = true;
    } else if (result.trans_count >= thresholds.min_cis_reads &&
               result.TransFraction() >= thresholds.min_cis_fraction) {
      result.is_phased = true;
      result.is_cis = false;
    }
  }

  return result;
}

namespace {

// Returns true if the pair is an insertion+SNV combination
bool IsInsertionSnvPair(const Variant& a, const Variant& b) {
  return (a.IsInsertion() && b.IsSnv()) || (a.IsSnv() && b.IsInsertion());
}

// A candidate variant pair (indices into the variant list) to phase.
struct PairTask {
  s32 vi;
  s32 vj;
};

// Half-open range [start, end) of task indices assigned to one worker.
struct ChunkRange {
  size_t start;
  size_t end;
};

// Enumerate all candidate variant pairs across groups, applying the distance
// and insertion+SNV filters. Returns the tasks; logs any skipped pairs.
vec<PairTask> EnumeratePairTasks(const vec<Variant>& variants, const vec<PhaseGroup>& groups, const Params& params) {
  vec<PairTask> tasks;
  s32 skipped_ins_snv = 0;
  for (const auto& group : groups) {
    for (size_t i = 0; i < group.variant_indices.size(); ++i) {
      for (size_t j = i + 1; j < group.variant_indices.size(); ++j) {
        const s32 vi = group.variant_indices[i];
        const s32 vj = group.variant_indices[j];
        const auto uvi = static_cast<size_t>(vi);
        const auto uvj = static_cast<size_t>(vj);
        if (std::abs(variants[uvi].pos - variants[uvj].pos) > params.max_phase_distance) {
          continue;
        }
        // Skip indel+SNV pairs when requested (high FP rate near insertions)
        if (params.skip_ins_snv && IsInsertionSnvPair(variants[uvi], variants[uvj])) {
          ++skipped_ins_snv;
          continue;
        }
        tasks.push_back({vi, vj});
      }
    }
  }
  if (skipped_ins_snv > 0) {
    Logging::Info("Skipped {} insertion+SNV pairs (--skip-ins-snv)", skipped_ins_snv);
  }
  return tasks;
}

// Warn if the requested tumor read group matches no @RG line in the BAM header.
// Reads are selected by an exact RG match or a "<tumor_rg>-" prefix match
// (Mutect2 bamout suffixes). When nothing matches, every read is skipped and no
// pair can be phased, so surface the likely cause up front.
void WarnIfTumorRgAbsent(const io::AlignmentReader& bam_handle, const std::string& tumor_rg) {
  if (bam_handle.hdr == nullptr || tumor_rg.empty()) {
    return;
  }
  sam_hdr_t* const hdr = bam_handle.hdr.get();
  const s32 n_rg = sam_hdr_count_lines(hdr, "RG");
  if (n_rg <= 0) {
    Logging::Warn("BAM header has no @RG lines; no reads will match --tumor-sample-name '{}'", tumor_rg);
    return;
  }

  const std::string tumor_rg_prefix = tumor_rg + "-";
  vec<std::string> rg_ids;
  bool matched = false;
  for (s32 i = 0; i < n_rg; ++i) {
    kstring_t id_val = KS_INITIALIZE;
    if (sam_hdr_find_tag_pos(hdr, "RG", i, "ID", &id_val) == 0 && id_val.s != nullptr) {
      const std::string_view rg_id(id_val.s);
      rg_ids.emplace_back(rg_id);
      if (rg_id == tumor_rg || rg_id.starts_with(tumor_rg_prefix)) {
        matched = true;
      }
    }
    ks_free(&id_val);
  }

  if (!matched) {
    Logging::Warn(
        "--tumor-sample-name '{}' matches none of the {} read group(s) in the BAM ([{}]); "
        "all reads will be skipped and no variants can be phased. "
        "Pass the read group ID (@RG ID), not the sample name (SM).",
        tumor_rg,
        n_rg,
        fmt::join(rg_ids, ", "));
  }
}

// Phase the pair tasks in [range.start, range.end) using a thread-local BAM handle,
// writing each result into its own slot in `results` (no locking needed).
void PhasePairChunk(const fs::path& bam_path,
                    const std::string& tumor_rg,
                    const vec<Variant>& variants,
                    const Params& params,
                    const vec<PairTask>& tasks,
                    const ChunkRange range,
                    vec<PhasingResult>& results) {
  const io::AlignmentReader bam_handle = io::OpenAlignmentReader(bam_path);
  if (bam_handle.bam_file == nullptr || bam_handle.hdr == nullptr || bam_handle.idx == nullptr) {
    Logging::Error("Thread failed to open BAM: {}", bam_path.string());
    return;
  }
  for (size_t i = range.start; i < range.end; ++i) {
    const auto& task = tasks[i];
    results[i] = PhaseVariantPair(bam_handle,
                                  tumor_rg,
                                  variants[static_cast<size_t>(task.vi)],
                                  variants[static_cast<size_t>(task.vj)],
                                  {.min_shared_reads = params.min_shared_reads,
                                   .min_cis_reads = params.min_cis_reads,
                                   .min_cis_fraction = params.min_cis_fraction,
                                   .min_base_quality = params.min_base_quality});
    results[i].variant_idx_a = task.vi;
    results[i].variant_idx_b = task.vj;
  }
}

}  // namespace

vec<PhaseBlock> RunPhasing(const fs::path& bam_path,
                           const std::string& tumor_rg,
                           const vec<Variant>& variants,
                           const vec<PhaseGroup>& groups,
                           const Params& params) {
  vec<PhasingResult> all_results;

  const io::AlignmentReader bam_handle = io::OpenAlignmentReader(bam_path);
  if (bam_handle.bam_file == nullptr || bam_handle.hdr == nullptr || bam_handle.idx == nullptr) {
    Logging::Error("Failed to open BAM: {}", bam_path.string());
    return {};
  }
  WarnIfTumorRgAbsent(bam_handle, tumor_rg);

  s32 skipped_ins_snv = 0;
  for (const auto& group : groups) {
    // Phase all pairs within the group
    for (size_t i = 0; i < group.variant_indices.size(); ++i) {
      for (size_t j = i + 1; j < group.variant_indices.size(); ++j) {
        const s32 vi = group.variant_indices[i];
        const s32 vj = group.variant_indices[j];
        const auto uvi = static_cast<size_t>(vi);
        const auto uvj = static_cast<size_t>(vj);

        // Skip pairs beyond max distance
        if (std::abs(variants[uvi].pos - variants[uvj].pos) > params.max_phase_distance) {
          continue;
        }

        // Skip indel+SNV pairs when requested (high FP rate near insertions)
        if (params.skip_ins_snv && IsInsertionSnvPair(variants[uvi], variants[uvj])) {
          ++skipped_ins_snv;
          continue;
        }

        auto result = PhaseVariantPair(bam_handle,
                                       tumor_rg,
                                       variants[uvi],
                                       variants[uvj],
                                       {.min_shared_reads = params.min_shared_reads,
                                        .min_cis_reads = params.min_cis_reads,
                                        .min_cis_fraction = params.min_cis_fraction,
                                        .min_base_quality = params.min_base_quality});
        result.variant_idx_a = vi;
        result.variant_idx_b = vj;

        if (result.is_phased) {
          Logging::Debug("Phased: {}:{} <-> {}:{} cis={} trans={} ref={}",
                         variants[uvi].chrom,
                         variants[uvi].pos + 1,
                         variants[uvj].chrom,
                         variants[uvj].pos + 1,
                         result.cis_count,
                         result.trans_count,
                         result.ref_count);
        }

        all_results.push_back(result);
      }
    }
  }

  auto blocks = BuildPhaseBlocks(all_results, variants);

  if (skipped_ins_snv > 0) {
    Logging::Info("Skipped {} insertion+SNV pairs (--skip-ins-snv)", skipped_ins_snv);
  }
  Logging::Info("Found {} phase blocks from {} phased pairs",
                blocks.size(),
                std::ranges::count_if(all_results, [](const PhasingResult& r) { return r.is_phased; }));

  return blocks;
}

vec<PhaseBlock> RunPhasingParallel(const fs::path& bam_path,
                                   const std::string& tumor_rg,
                                   const vec<Variant>& variants,
                                   const vec<PhaseGroup>& groups,
                                   const Params& params) {
  if (params.threads <= 1 || groups.size() <= 1) {
    return RunPhasing(bam_path, tumor_rg, variants, groups, params);
  }

  // Validate the tumor read group once before spawning threads.
  {
    const io::AlignmentReader bam_handle = io::OpenAlignmentReader(bam_path);
    if (bam_handle.bam_file == nullptr || bam_handle.hdr == nullptr || bam_handle.idx == nullptr) {
      Logging::Error("Failed to open BAM: {}", bam_path.string());
      return {};
    }
    WarnIfTumorRgAbsent(bam_handle, tumor_rg);
  }

  // Enumerate all candidate pairs for fine-grained parallelism
  const vec<PairTask> tasks = EnumeratePairTasks(variants, groups, params);
  Logging::Info("Parallelizing {} variant pairs across {} threads", tasks.size(), params.threads);

  // Each task gets its own result slot — no locking needed
  vec<PhasingResult> results(tasks.size());

  tf::Executor executor(params.threads);
  tf::Taskflow taskflow;

  // Partition pairs into chunks for each thread.
  // Each thread opens its own BAM handle (htslib is not thread-safe per handle).
  const size_t n_threads = params.threads;
  const size_t chunk_size = (tasks.size() + n_threads - 1) / n_threads;

  for (size_t t = 0; t < n_threads; ++t) {
    const size_t start = t * chunk_size;
    const size_t end = std::min(start + chunk_size, tasks.size());
    if (start >= tasks.size()) {
      break;
    }
    taskflow.emplace([&bam_path, &tumor_rg, &variants, &params, &tasks, &results, start, end]() {
      PhasePairChunk(bam_path, tumor_rg, variants, params, tasks, {.start = start, .end = end}, results);
    });
  }

  executor.run(taskflow).wait();

  // Merge results and build phase blocks with haplotype assignment
  const vec<PhaseBlock> blocks = BuildPhaseBlocks(results, variants);

  const auto phased_count = std::ranges::count_if(results, [](const PhasingResult& r) { return r.is_phased; });
  Logging::Info("Found {} phase blocks from {} phased pairs ({} threads)", blocks.size(), phased_count, params.threads);

  if (phased_count == 0 && !results.empty()) {
    const auto pairs_with_reads =
        std::ranges::count_if(results, [](const PhasingResult& r) { return r.total_spanning > 0; });
    if (pairs_with_reads == 0) {
      Logging::Warn(
          "No candidate pair had any spanning reads. Check that --bam-input contains the reads "
          "(e.g. the Mutect2 bamout) and that --tumor-sample-name matches an @RG ID in the BAM.");
    } else {
      Logging::Info("{} of {} candidate pairs had spanning reads but insufficient concordant support to phase",
                    pairs_with_reads,
                    results.size());
    }
  }

  return blocks;
}

}  // namespace xoos::svc::phaser
