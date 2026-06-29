#include "duplex-match.h"

#include <xoos/error/error.h>
#include <xoos/log/logging.h>

#include <algorithm>

#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

static void UpdateIndices(const std::vector<gtl::flat_hash_map<u64, BarcodeMatch>>& hash_table, size_t& min_index,
                          size_t& max_index) {
  for (size_t i = 0; i < hash_table.size(); ++i) {
    const auto& entry = hash_table[i];
    if (!entry.empty()) {
      if (min_index > i) {
        min_index = i;
      }
      if (max_index < i) {
        max_index = i;
      }
    }
  }
}

/**
 * @brief Based on other constructor written but generic and only has 2 SeqMatcher objects.
 * @param seqs_5p SeqMatcher object for 5' sequences.
 * @param seqs_3p SeqMatcher object for 3' sequences.
 * @param barcode_type_5p Barcode type for 5' sequences.
 * @param barcode_type_3p Barcode type for 3' sequences.
 */
CascadedLUTs::CascadedLUTs(const SeqMatcher& seqs_5p,  // NOLINT - we initialize all fields in the constructor
                           const SeqMatcher& seqs_3p, DuplexMatch::BarcodeType barcode_type_5p,
                           DuplexMatch::BarcodeType barcode_type_3p) {
  const auto& hash_tables_5p{seqs_5p.Lut()->HashTables()};
  const auto& hash_tables_3p{seqs_3p.Lut()->HashTables()};

  u64 max_5p = 0;
  u64 min_5p = hash_tables_5p.size();
  UpdateIndices(hash_tables_5p, min_5p, max_5p);

  u64 max_3p = 0;
  u64 min_3p = hash_tables_3p.size();
  UpdateIndices(hash_tables_3p, min_3p, max_3p);

  // The minimum length of all hash tables determines the number of entries in the first LUT. For every
  // possible bit combination (kLut0Bits), we need to reserve multiple entries to reflect the possible
  // barcode types.
  _lut0.resize((1u << kLut0Bits) * static_cast<std::size_t>(DuplexMatch::BarcodeType::kUnknown), 0);

  // Allocate space for the "zero" entries.
  size_t max_length = std::max({max_5p, max_3p});
  size_t all_zero_size{1u << (2 * max_length - kLut0Bits)};
  AddMemory(0, all_zero_size);
  for (size_t i = 1; i < DuplexMatch::BarcodeType::kUnknown; ++i) {
    _lut0[i] = _lut0[0];  // Reuse the memory for the dummy entries
  }

  AddToLUTs(seqs_5p, barcode_type_5p);
  AddToLUTs(seqs_3p, barcode_type_3p);
  Logging::Debug("Size of LUT entries: {} bytes", MemoryUsage());
}

void CascadedLUTs::AddMemory(
    size_t index,
    size_t nr_entries) {  // This function adds memory to the second LUT where the actual LUT information is stored.
  _lut0[index] = static_cast<u32>(_lut1.size());
  _lut1.resize(_lut1.size() + nr_entries);
}

static std::string HashString(u64 hashvalue, u64 length) {
  std::string hash_str(length, '\0');
  for (size_t i = 0; i < length; ++i) {
    switch (hashvalue & kLast2BitsMask) {
      case k00:
        hash_str[i] = 'A';
        break;
      case k01:
        hash_str[i] = 'C';
        break;
      case k10:
        hash_str[i] = 'G';
        break;
      default:
        hash_str[i] = 'T';
        break;
    }
    hashvalue >>= kRightShift2;
  }
  return hash_str;
}

void CascadedLUTs::AddToLUTs(const SeqMatcher& seqs, DuplexMatch::BarcodeType barcode_type) {
  const auto& hash_tables = seqs.Lut()->HashTables();

  // (re)calculate nominal barcode length
  size_t min_length = hash_tables.size(), max_length = 0;
  for (size_t length = 0; length < hash_tables.size(); ++length) {
    const auto& hash_table = hash_tables[length];
    if (!hash_table.empty()) {
      if (min_length > length) {
        min_length = length;
      }
      if (max_length < length) {
        max_length = length;
      }
    }
  }

  // calculate nominal length
  const auto nominal_length = static_cast<s32>((min_length + max_length) >> 1);
  _nominal_length[barcode_type] = nominal_length;
  _max_length[barcode_type] = max_length;

  for (size_t length = 0; length < hash_tables.size(); ++length) {
    const auto& hash_table = hash_tables[length];
    if (!hash_table.empty()) {
      u16 l{DuplexMatch::kNone};
      switch (static_cast<int>(length) - nominal_length) {
        case -2:
          l = DuplexMatch::kMinus2;
          break;
        case -1:
          l = DuplexMatch::kMinus1;
          break;
        case 0:
          l = DuplexMatch::kZero;
          break;
        case 1:
          l = DuplexMatch::kPlus1;
          break;
        case 2:
          l = DuplexMatch::kPlus2;
          break;
        default:
          break;
      }

      DuplexMatch lut_entry;
      lut_entry.length = l;
      lut_entry.barcode_type = barcode_type;

      const auto lut1_size{1u << 2 * (max_length - kNrBasesLUT0)};

      for (const auto& [hash_value, barcode_match] : hash_table) {
        // Calculate the number of bases stored in the LUT. Note that the number of bases stored in the first
        // LUT might be larger than the minimum number of bases of the marker; in that case, we'll have to
        // add dummy entries.
        // Let's assume 11 bases in the first LUT.
        // SID nominal length is 12, let's consider distance 2 so minimum size is 10.
        // SID length 10: LUT0: BBBBBBBBBB* (10 bases + 1 wildcard), LUT1: *** (3 wildcards)
        // SID length 11: LUT0: BBBBBBBBBBB (11 bases), LUT1: *** (3 wildcards)
        // SID length 12: LUT0: BBBBBBBBBBB (11 bases), LUT1: B** (1 base, 2 wildcards)
        // SID length 13: LUT0: BBBBBBBBBBB (11 bases), LUT1: BB* (2 bases, 1 wildcard)
        // SID length 14: LUT0: BBBBBBBBBBB (11 bases), LUT1: BBB (3 bases)
        lut_entry.edist = barcode_match.edist;

        const size_t nr_wildcard_bases{max_length - length};
        const size_t nr_wildcard_bits{nr_wildcard_bases << 1};
        const size_t total_entries_wildcards{1u << nr_wildcard_bits};
        for (size_t i = 0; i < total_entries_wildcards; ++i) {
          // calculate the hash value for the wildcard entry
          size_t hash_value_wildcard = hash_value | (i << (2 * length));
          size_t index0{hash_value_wildcard & kLut0Mask};
          size_t index1{hash_value_wildcard >> kLut0Bits};
          const size_t index2{index0 * static_cast<std::size_t>(DuplexMatch::BarcodeType::kUnknown) +
                              static_cast<std::size_t>(barcode_type)};
          if (_lut0[index2] == 0) {
            AddMemory(index2, lut1_size);
          }

          auto& current_entry = _lut1[_lut0[index2] + index1];

          // Error checking: if the entry is already filled by a different SID or barcode type, then we have a problem.
          bool current_is_sid =
              current_entry.barcode_type == DuplexMatch::kSID5p || current_entry.barcode_type == DuplexMatch::kSID3p;
          bool new_is_sid = barcode_type == DuplexMatch::kSID5p || barcode_type == DuplexMatch::kSID3p;
          if (current_is_sid && new_is_sid && current_entry.barcode_id != barcode_match.barcode_id &&
              current_entry.length == lut_entry.length) {
            // Prefixing bit field with a unaray plus as suggested here: https://github.com/fmtlib/fmt/issues/1284
            throw error::Error("Input {} matches both SID {} and {} with edit distance {} and {} respectively.",
                               HashString(hash_value, length), +current_entry.barcode_id, +barcode_match.barcode_id,
                               +barcode_match.edist, +current_entry.edist);
          }

          if (lut_entry < current_entry) {
            // If the lut entry is a better match compared to current entry (smaller edit distance, longer sequence),
            // then replace current entry. Otherwise, keep the current entry.
            current_entry = lut_entry;
            current_entry.barcode_id = barcode_match.barcode_id;
          }
        }
      }
    }
  }
}

}  // namespace xoos::demux
