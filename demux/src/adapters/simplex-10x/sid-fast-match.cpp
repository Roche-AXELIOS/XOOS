#include "sid-fast-match.h"

#include <algorithm>
#include <array>
#include <cstring>

#include <xoos/log/logging.h>
#include <xoos/util/sequence-functions.h>

// SidFastMatch — hash-table-based fast-path for SID barcode assignment.
//
// Most SIMPLEX-10X reads (~86%) start with an exact runway at offset 0 or 1.
// For these reads we skip the full LUT + DFA + Bitap slow path and resolve
// the SID with a single hash-table lookup:
//
//   1. Match the runway (truncated or full) at the read start.
//   2. Encode the next 16 bases (SID + stem prefix) as a u32 key.
//   3. Look up the key — if unambiguous, return the SID assignment.
//
// The hash table is built by Build(): exact keys (edist=0) plus all
// 1-substitution variants (edist=1), with collisions marked ambiguous.
// Reads that miss or hit an ambiguous entry fall through to the slow path.

namespace xoos::demux {

// --- Encoding ---

std::optional<u32> SidFastMatch::Encode(const char* const seq,  // NOSONAR(S6188) -- raw pointer for hot-path arithmetic
                                        const u32 len) {
  u32 val = 0;
  for (u32 i = 0; i < len; ++i) {
    const auto bits =
        sequence::kBaseToBin[static_cast<u8>(seq[i])];  // NOSONAR(S810,M23_058) -- matches kBaseToBin convention
    if (bits > kMaxValidBase) {
      return std::nullopt;
    }
    val |= (static_cast<u32>(bits) << (i * kBitsPerBase));
  }
  return val;
}

// --- Construction ---

SidFastMatch::SidFastMatch(const std::string_view runway,
                           const u32 sid_len,
                           const u32 key_len,
                           gtl::flat_hash_map<u32, SidEntry> sid_map)
    : _runway(runway), _sid_len(sid_len), _key_len(key_len), _sid_map(std::move(sid_map)) {
}

std::optional<u32> SidFastMatch::ValidatePool(const BarcodePool& sid_pool) {
  const auto sid_len = static_cast<u32>(sid_pool.front().sequence.size());

  for (const auto& sid : sid_pool) {
    if (sid.sequence.size() != sid_len) {
      Logging::Warn("SidFastMatch: variable SID lengths ({} vs {}), skipping fast-path", sid_len, sid.sequence.size());
      return std::nullopt;
    }
  }

  if (sid_len > kMaxEncodableLen) {
    Logging::Warn(
        "SidFastMatch: SID length {} exceeds max encodable length {}, skipping fast-path", sid_len, kMaxEncodableLen);
    return std::nullopt;
  }

  for (const auto& sid : sid_pool) {
    if (!Encode(sid.sequence.data(), sid_len)) {
      Logging::Warn("SidFastMatch: non-ACGT base in SID '{}', skipping fast-path", sid.sequence);
      return std::nullopt;
    }
  }

  return sid_len;
}

void SidFastMatch::InsertMismatchVariants(const BarcodePool& sid_pool,
                                          const u32 sid_len,
                                          const u32 stem_bits,
                                          gtl::flat_hash_map<u32, SidEntry>& sid_map) {
  // Mutate SID positions 0..(sid_len-2), excluding the last SID base and
  // the stem prefix. The last SID base borders the stem prefix — a
  // sequencing error there is indistinguishable from a 1bp boundary shift,
  // which affects the trim position. The fast path cannot detect this, so
  // these reads are deferred to the slow path for full scoring.
  //
  // Collision handling:
  //   Case 1: variant of SID_A collides with exact match of SID_B (edist=0).
  //     The exact match takes priority; do not overwrite.
  //   Case 2: variant of SID_A collides with variant of SID_B (edist=1).
  //     Assignment is ambiguous — mark as kAmbiguous (falls through to slow path).
  constexpr u32 kBaseCount = 4;
  constexpr u32 kBaseMask = 0x3;
  const u32 mismatch_end = sid_len > 0 ? sid_len - 1 : 0;

  for (const auto& sid : sid_pool) {
    const u32 full_key = *Encode(sid.sequence.data(), sid_len) | stem_bits;
    for (u32 pos = 0; pos < mismatch_end; ++pos) {
      const u32 shift = pos * kBitsPerBase;
      const u32 orig_bits = (full_key >> shift) & kBaseMask;
      for (u32 alt = 0; alt < kBaseCount; ++alt) {
        if (alt == orig_bits) {
          continue;
        }
        const u32 variant_key = (full_key & ~(kBaseMask << shift)) | (alt << shift);
        const auto [it, inserted] = sid_map.emplace(variant_key, SidEntry{sid.id, 1});
        if (!inserted && it->second.sid_id != sid.id && it->second.edist != 0) {
          // Two cases when `inserted` is false (key already exists):
          //  1. Different SID, existing is exact match (edist == 0) → keep
          //     exact match; it was inserted in phase 1 and takes priority.
          //  2. Different SID, both at edist == 1 → ambiguous, mark below.
          it->second.sid_id = kAmbiguous;
        }
      }
    }
  }
}

std::optional<SidFastMatch> SidFastMatch::Build(const std::string_view runway,
                                                const BarcodePool& sid_pool,
                                                const std::string_view stem) {
  u32 stem_prefix_len = kStemPrefixLen;
  // The truncated runway is runway[1:], so we need at least 2 characters.
  constexpr u32 kMinRunwayLen = 2;
  if (runway.size() < kMinRunwayLen) {
    Logging::Warn("SidFastMatch: runway too short (need >= {} for truncated form), skipping fast-path", kMinRunwayLen);
    return std::nullopt;
  }
  if (sid_pool.empty()) {
    Logging::Warn("SidFastMatch: empty SID pool, skipping fast-path");
    return std::nullopt;
  }

  const auto validated_sid_len = ValidatePool(sid_pool);
  if (!validated_sid_len) {
    return std::nullopt;
  }
  const u32 sid_len = *validated_sid_len;

  // --- Stem prefix adjustment (non-critical, silently truncated) ---

  stem_prefix_len = std::min(stem_prefix_len, static_cast<u32>(stem.size()));

  const u32 max_prefix = kMaxEncodableLen - sid_len;
  if (stem_prefix_len > max_prefix) {
    Logging::Info(
        "SidFastMatch: truncating stem prefix from {} to {} (SID length {})", stem_prefix_len, max_prefix, sid_len);
    stem_prefix_len = max_prefix;
  }

  u32 stem_bits = 0;
  if (stem_prefix_len > 0) {
    auto encoded = Encode(stem.data(), stem_prefix_len);
    if (!encoded) {
      Logging::Info("SidFastMatch: non-ACGT base in stem prefix, building without stem prefix");
      stem_prefix_len = 0;
    } else {
      stem_bits = *encoded << (sid_len * kBitsPerBase);
    }
  }

  const u32 effective_key_len = sid_len + stem_prefix_len;

  // --- Hash-table population ---

  gtl::flat_hash_map<u32, SidEntry> sid_map;

  // Phase 1: insert exact matches (edist=0).
  for (const auto& sid : sid_pool) {
    const u32 full_key = *Encode(sid.sequence.data(), sid_len) | stem_bits;
    const auto [it, inserted] = sid_map.emplace(full_key, SidEntry{sid.id, 0});
    if (!inserted && it->second.sid_id != sid.id) {
      Logging::Warn("SidFastMatch: exact-match key collision between SID {} and SID {}, skipping fast-path",
                    it->second.sid_id,
                    sid.id);
      return std::nullopt;
    }
  }

  // Phase 2: insert 1-substitution variants (edist=1).
  InsertMismatchVariants(sid_pool, sid_len, stem_bits, sid_map);

  Logging::Info("SidFastMatch: built fast-path with {} SIDs, key_len={} (SID {} + stem prefix {}), {} hash entries",
                sid_pool.size(),
                effective_key_len,
                sid_len,
                stem_prefix_len,
                sid_map.size());

  return SidFastMatch(runway, sid_len, effective_key_len, std::move(sid_map));
}

// --- Lookup ---

std::optional<SidFastMatch::Result> SidFastMatch::TryMatch(const char* const seq, const size_t length) const {
  const auto full_len = static_cast<u32>(_runway.size());
  const auto trunc_len = full_len - 1;

  // Try truncated runway first — it's the more common form (~50% of reads
  // lose the leading base during sequencing vs ~40% that retain it).
  struct Attempt {
    const char* runway_ptr;
    u32 runway_len;
  };

  // Truncated runway (e.g. AACAA) first — more common (~50% of reads lose
  // the leading base). Full runway (e.g. CAACAA) second (~40%).
  const std::array<Attempt, 2> attempts = {{
      {_runway.data() + 1, trunc_len},
      {_runway.data(), full_len},
  }};

  for (const auto& [runway_ptr, runway_len] : attempts) {
    const auto min_len = runway_len + _key_len;
    if (min_len > length) {
      continue;
    }

    // Step 1: exact runway match.
    if (std::memcmp(seq, runway_ptr, runway_len) != 0) {  // NOSONAR(S5356) -- char* to void* is safe for memcmp
      continue;
    }

    // Step 2: encode SID + stem prefix as a single u32 key.
    const u32 key_start = runway_len;
    const auto encoded = Encode(seq + key_start, _key_len);
    if (!encoded) {
      // Non-ACGT base in the key region.
      continue;
    }

    // Step 3: hash-table lookup.
    const auto it = _sid_map.find(*encoded);
    if (it == _sid_map.end()) {
      // No matching SID within 1 substitution.
      continue;
    }
    if (it->second.sid_id == kAmbiguous) {
      // Ambiguous 1-mismatch collision — defer to slow path.
      continue;
    }

    // Step 4: return the SID assignment.
    // insert_start points right after the SID (stem remains in the insert).
    return Result{it->second.sid_id, it->second.edist, key_start + _sid_len};
  }

  return std::nullopt;
}

}  // namespace xoos::demux
