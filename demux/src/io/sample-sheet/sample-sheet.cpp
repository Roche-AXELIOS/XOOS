#include "sample-sheet.h"

#include <algorithm>
#include <fstream>
#include <locale>
#include <optional>
#include <sstream>
#include <string_view>

#include <fmt/ranges.h>

#include <csv.hpp>

#include <xoos/error/error.h>
#include <xoos/log/logging.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>
#include <xoos/util/sequence-functions.h>

#include "adapter-design/adapter-design.h"
#include "task/flow-context.h"

namespace xoos::demux {

namespace {

vec<std::string> SplitHeaderLine(const std::string& line) {
  const char delimiter = line.contains('\t') ? '\t' : ',';
  vec<std::string> result;
  std::stringstream stream(line);
  while (true) {
    std::string field;
    if (!std::getline(stream, field, delimiter)) {
      break;
    }
    if (!field.empty() && field.back() == '\r') {
      field.pop_back();
    }
    result.emplace_back(std::move(field));
  }
  return result;
}

bool ContainsHeader(const vec<std::string>& headers, const std::string_view name) {
  return std::ranges::contains(headers, name);
}

void ValidateHeader(const fs::path& path) {
  std::ifstream file(path);
  std::string line;
  if (!std::getline(file, line)) {
    throw error::Error("Sample sheet '{}' is empty", path.string());
  }

  const auto headers = SplitHeaderLine(line);
  if (ContainsHeader(headers, kSampleSidColumnHeader) && ContainsHeader(headers, kSampleNameColumnHeader)) {
    return;
  }

  throw error::Error("Sample sheet '{}' must contain '{}' and '{}' columns",
                     path.string(),
                     kSampleNameColumnHeader,
                     kSampleSidColumnHeader);
}

std::string NormalizeSidSequence(std::string sid_sequence, const std::string_view field_name) {
  const auto& locale = std::locale::classic();
  std::ranges::transform(
      sid_sequence, sid_sequence.begin(), [&locale](const char base) { return std::toupper(base, locale); });
  for (const char base : sid_sequence) {
    if (!sequence::IsACGT(base)) {
      throw error::Error("Sample sheet field '{}' contains invalid SID character '{}'", field_name, base);
    }
  }
  return sid_sequence;
}

std::pair<std::string, std::string> ParseSampleSid(std::string sample_sid) {
  const auto separator = sample_sid.find('-');
  if (separator == std::string::npos) {
    sample_sid = NormalizeSidSequence(std::move(sample_sid), kSampleSidColumnHeader);
    return {sample_sid, sample_sid};
  }

  if (separator == 0 || separator == sample_sid.size() - 1 ||
      sample_sid.find('-', separator + 1) != std::string::npos) {
    throw error::Error("Sample sheet field '{}' must be either SID or SID5-SID3", kSampleSidColumnHeader);
  }

  auto sid_5p = NormalizeSidSequence(sample_sid.substr(0, separator), kSampleSidColumnHeader);
  auto sid_3p = NormalizeSidSequence(sample_sid.substr(separator + 1), kSampleSidColumnHeader);
  return {std::move(sid_5p), std::move(sid_3p)};
}

std::string FormatSidPair(const SampleSheetEntry& entry) {
  return entry.sid_5p + "-" + entry.sid_3p;
}

std::string FormatSidForMessage(const SampleSheetEntry& entry) {
  if (entry.sid_5p == entry.sid_3p) {
    return entry.sid_5p;
  }
  return FormatSidPair(entry);
}

DualSidPartialMatch ApplyDualSidMode(DualSidPartialMatch match, const DualSidMode mode) {
  using enum DualSidMode;
  switch (mode) {
    case kNone:
      break;
    case kNeed5p:
      match.allow_3p_partial_match = false;
      break;
    case kNeed3p:
      match.allow_5p_partial_match = false;
      break;
    case kNeedBoth:
      match.allow_5p_partial_match = false;
      match.allow_3p_partial_match = false;
      break;
    default:
      break;
  }
  return match;
}

u32 FindOrAddSideSid(BarcodePool& pool, StrUnorderedMap<u32>& id_by_sequence, const std::string& sequence) {
  const auto it = id_by_sequence.find(sequence);
  if (it != id_by_sequence.end()) {
    return it->second;
  }

  const auto id = static_cast<u32>(pool.size());
  id_by_sequence[sequence] = id;
  pool.emplace_back(id, sequence, sequence);
  return id;
}

char AsciiToUpper(const char base) {
  return std::toupper(base, std::locale::classic());
}

std::string ApplyTransform(std::string sequence, const SeqTransform transform) {
  switch (transform) {
    using enum SeqTransform;
    case kReverse:
      std::ranges::reverse(sequence);
      return sequence;
    case kReverseComplement:
      return sequence::ReverseComplement(sequence);
    case kNone:
      return sequence;
    default:
      throw error::Error("Unsupported SID sequence transform");
  }
}

struct LegacySidSignature {
  SeqTransform transform{SeqTransform::kNone};
  std::string read_prefix;
  std::string read_suffix;

  friend bool operator==(const LegacySidSignature&, const LegacySidSignature&) = default;
};

struct LegacySidMatch {
  std::string canonical;
  LegacySidSignature signature;

  friend bool operator==(const LegacySidMatch&, const LegacySidMatch&) = default;
};

bool MatchesStructuralFlanks(const std::string_view supplied_prefix,
                             const std::string_view supplied_suffix,
                             const MatcherComponentDefinition& matcher) {
  return supplied_prefix.size() <= matcher.structural_prefix.size() &&
         std::string_view{matcher.structural_prefix}.ends_with(supplied_prefix) &&
         supplied_suffix.size() <= matcher.structural_suffix.size() &&
         std::string_view{matcher.structural_suffix}.starts_with(supplied_suffix);
}

std::optional<LegacySidMatch> FindLegacySidMatch(const std::string_view supplied,
                                                 const BarcodePool& canonical_sids,
                                                 const AdapterDesign& design) {
  vec<LegacySidMatch> candidates;
  for (const AdapterComponentDefinition& component : design.components) {
    const auto* const matcher = std::get_if<MatcherComponentDefinition>(&component);
    if (matcher == nullptr || matcher->kind != MatcherKind::kSid ||
        (matcher->structural_prefix.empty() && matcher->structural_suffix.empty())) {
      continue;
    }

    for (const Barcode& barcode : canonical_sids) {
      const std::string oriented = ApplyTransform(barcode.sequence, matcher->transform);
      if (supplied.size() <= oriented.size()) {
        continue;
      }

      std::size_t position = supplied.find(oriented);
      while (position != std::string::npos) {
        const std::string_view prefix{supplied.data(), position};
        const std::size_t suffix_position = position + oriented.size();
        const std::string_view suffix{supplied.data() + suffix_position, supplied.size() - suffix_position};
        if ((!prefix.empty() || !suffix.empty()) && MatchesStructuralFlanks(prefix, suffix, *matcher)) {
          const LegacySidMatch candidate{
              barcode.sequence,
              {matcher->transform, std::string{prefix}, std::string{suffix}},
          };
          if (std::ranges::find(candidates, candidate) == candidates.end()) {
            candidates.push_back(candidate);
          }
        }
        position = supplied.find(oriented, position + 1);
      }
    }
  }
  return candidates.size() == 1U ? std::optional<LegacySidMatch>{candidates.front()} : std::nullopt;
}

struct SidPoolAnalysis {
  size_t canonical_count{0};
  size_t reverse_complement_count{0};
  vec<std::string> legacy_canonical_sequences;
  std::optional<LegacySidSignature> legacy_signature;
  bool legacy_signature_is_uniform{true};
  vec<std::string> unrecognized_sid_sample_names;
};

void RecordLegacyMatch(SidPoolAnalysis& analysis, const LegacySidMatch& match) {
  if (!analysis.legacy_signature.has_value()) {
    analysis.legacy_signature = match.signature;
  } else if (analysis.legacy_signature.value() != match.signature) {
    analysis.legacy_signature_is_uniform = false;
  }
  analysis.legacy_canonical_sequences.push_back(match.canonical);
}

SidPoolAnalysis AnalyzeSidPool(const SidPool& sample_sids,
                               const BarcodePool& canonical_sids,
                               const StrUnorderedSet& canonical,
                               const AdapterDesign& design) {
  SidPoolAnalysis analysis;
  for (const Barcode& barcode : sample_sids) {
    if (canonical.contains(barcode.sequence)) {
      ++analysis.canonical_count;
    } else if (canonical.contains(sequence::ReverseComplement(barcode.sequence))) {
      ++analysis.reverse_complement_count;
    } else if (const auto legacy_match = FindLegacySidMatch(barcode.sequence, canonical_sids, design);
               legacy_match.has_value()) {
      RecordLegacyMatch(analysis, *legacy_match);
    } else {
      analysis.unrecognized_sid_sample_names.push_back(barcode.name);
    }
  }
  return analysis;
}

bool HasMixedRepresentations(const SidPoolAnalysis& analysis) {
  const size_t legacy_count = analysis.legacy_canonical_sequences.size();
  const bool has_canonical = analysis.canonical_count > 0U;
  const bool has_reverse_complement = analysis.reverse_complement_count > 0U;
  const bool has_legacy = legacy_count > 0U;
  const bool has_mixed_legacy_signature = has_legacy && !analysis.legacy_signature_is_uniform;
  if (has_mixed_legacy_signature) {
    return true;
  }
  if (has_canonical) {
    return has_reverse_complement || has_legacy;
  }
  return has_reverse_complement && has_legacy;
}

vec<std::string> NormalizeSidCaseInPlace(SidPool& sample_sids) {
  vec<std::string> converted_sample_names;
  for (Barcode& barcode : sample_sids) {
    std::string upper = barcode.sequence;
    std::ranges::transform(upper, upper.begin(), AsciiToUpper);
    if (upper != barcode.sequence) {
      converted_sample_names.push_back(barcode.name);
      barcode.sequence = std::move(upper);
    }
  }
  return converted_sample_names;
}

}  // namespace

void NormalizeSidPoolToCanonical(SidPool& sample_sids, const BarcodePool& canonical_sids, const AdapterDesign& design) {
  if (canonical_sids.empty() || sample_sids.empty()) {
    return;
  }

  const vec<std::string> converted_sample_names = NormalizeSidCaseInPlace(sample_sids);
  if (!converted_sample_names.empty()) {
    Logging::Warn("Normalized {} sample-sheet SID sequence(s) to uppercase before applying the adapter DSL: {}.",
                  converted_sample_names.size(),
                  fmt::join(converted_sample_names, ", "));
  }

  StrUnorderedSet canonical;
  canonical.reserve(canonical_sids.size());
  for (const Barcode& barcode : canonical_sids) {
    canonical.insert(barcode.sequence);
  }

  const SidPoolAnalysis analysis = AnalyzeSidPool(sample_sids, canonical_sids, canonical, design);
  const size_t legacy_count = analysis.legacy_canonical_sequences.size();
  if (analysis.reverse_complement_count == sample_sids.size()) {
    for (Barcode& barcode : sample_sids) {
      barcode.sequence = sequence::ReverseComplement(barcode.sequence);
    }
    Logging::Warn(
        "Legacy sample-sheet support: all {} SID(s) were supplied in reverse-complement orientation. "
        "Converting the entire sheet to canonical SID sequences before applying the adapter DSL. Consider updating "
        "the sample sheet to use canonical SID sequences.",
        sample_sids.size());
    return;
  }

  if (legacy_count == sample_sids.size() && analysis.legacy_signature_is_uniform) {
    for (std::size_t i = 0; i < sample_sids.size(); ++i) {
      sample_sids[i].sequence = analysis.legacy_canonical_sequences[i];
    }
    Logging::Warn(
        "Legacy sample-sheet support: all {} SID(s) use the same orientation with adapter structural bases. "
        "Converting the entire sheet to canonical SID sequences before applying adapter '{}'. Consider updating the "
        "sample sheet to use canonical SID sequences.",
        sample_sids.size(),
        design.name);
    return;
  }

  if (analysis.canonical_count == sample_sids.size()) {
    return;
  }
  if (HasMixedRepresentations(analysis)) {
    Logging::Warn(
        "Sample-sheet SIDs use mixed orientations or structural-flank representations; preserving them as supplied.");
  }
  if (!analysis.unrecognized_sid_sample_names.empty()) {
    Logging::Warn(
        "{} sample-sheet SID(s) are not recognized as canonical or supported legacy forms: {}. Preserving "
        "them as supplied.",
        analysis.unrecognized_sid_sample_names.size(),
        fmt::join(analysis.unrecognized_sid_sample_names, ", "));
  }
}

SampleSheetEntries ReadSampleSheet(const fs::path& path) {
  SampleSheetEntries entries;
  StrUnorderedMap<std::string> sid_pair2sample_name;
  StrUnorderedMap<std::string> sample_name2sid;
  ValidateHeader(path);
  csv::CSVFormat format;
  format.delimiter({',', '\t'});
  auto reader = csv::CSVReader(path.string(), format);
  for (const auto& row : reader) {
    auto sample_name = row[kSampleNameColumnHeader].get<std::string>();

    auto [sid_5p, sid_3p] = ParseSampleSid(row[kSampleSidColumnHeader].get<std::string>());
    auto entry = SampleSheetEntry{.sid_5p = std::move(sid_5p), .sid_3p = std::move(sid_3p), .name = sample_name};

    if (sample_name == kReservedRawFailedName) {
      throw error::Error("Sample name '{}' is reserved for unassigned reads and cannot be used in the sample sheet",
                         kReservedRawFailedName);
    }

    const auto sample_sid = FormatSidPair(entry);
    const auto it = sid_pair2sample_name.find(sample_sid);
    if (it != sid_pair2sample_name.end()) {
      throw error::Error("Found duplicate SID in sample sheet for '{}' and '{}'", sample_name, it->second);
    }
    sid_pair2sample_name[sample_sid] = sample_name;

    const auto it_name = sample_name2sid.find(sample_name);
    if (it_name != sample_name2sid.end()) {
      throw error::Error(
          "Found duplicate sample name in sample sheet for '{}' and '{}'", FormatSidForMessage(entry), it_name->second);
    }
    sample_name2sid[sample_name] = FormatSidForMessage(entry);

    entries.emplace_back(std::move(entry));
  }
  return entries;
}

SampleSheetSidPools LoadSampleSheetSidPools(const fs::path& sample_sheet, const DualSidMode dual_sid_mode) {
  const auto entries = ReadSampleSheet(sample_sheet);
  auto result = SampleSheetSidPools{};
  result.sid_pool.reserve(entries.size());
  result.matching_policy.partial_matches.reserve(entries.size());

  StrUnorderedMap<u32> sid_5p_counts;
  StrUnorderedMap<u32> sid_3p_counts;
  for (const auto& entry : entries) {
    ++sid_5p_counts[entry.sid_5p];
    ++sid_3p_counts[entry.sid_3p];
    result.has_dual_sids = result.has_dual_sids || entry.sid_5p != entry.sid_3p;
  }

  StrUnorderedMap<u32> sid_5p_id_by_sequence;
  StrUnorderedMap<u32> sid_3p_id_by_sequence;

  struct PairingWith3pId {
    u32 sid_3p_id;
    DualSidPairing pairing;
  };

  vec<PairingWith3pId> pairings;
  pairings.reserve(entries.size());
  u32 id = 0;
  for (const auto& entry : entries) {
    const auto sid_pair = FormatSidPair(entry);
    result.sid_pool.emplace_back(id, sid_pair, entry.name);
    const u32 sid_5p_id = FindOrAddSideSid(result.sid_5p_pool, sid_5p_id_by_sequence, entry.sid_5p);
    const u32 sid_3p_id = FindOrAddSideSid(result.sid_3p_pool, sid_3p_id_by_sequence, entry.sid_3p);
    if (sid_5p_counts.at(entry.sid_5p) == 1) {
      result.sid_5p_pool.at(sid_5p_id).name = entry.name;
    }
    if (sid_3p_counts.at(entry.sid_3p) == 1) {
      result.sid_3p_pool.at(sid_3p_id).name = entry.name;
    }
    result.matching_policy.sample_by_sid_pair[DualSidPairKey(sid_5p_id, sid_3p_id)] = id;
    pairings.emplace_back(PairingWith3pId{sid_3p_id, {.sid_5p_id = sid_5p_id, .sample_id = id}});

    result.matching_policy.partial_matches.emplace_back(
        ApplyDualSidMode(DualSidPartialMatch{.allow_5p_partial_match = sid_5p_counts.at(entry.sid_5p) == 1,
                                             .allow_3p_partial_match = sid_3p_counts.at(entry.sid_3p) == 1},
                         dual_sid_mode));
    ++id;
  }

  result.matching_policy.pairing_offsets_by_3p_sid.assign(result.sid_3p_pool.size() + 1, 0);
  for (const auto& pairing : pairings) {
    auto& pairing_count = result.matching_policy.pairing_offsets_by_3p_sid.at(pairing.sid_3p_id + 1);
    ++pairing_count;
  }
  for (std::size_t i = 1; i < result.matching_policy.pairing_offsets_by_3p_sid.size(); ++i) {
    result.matching_policy.pairing_offsets_by_3p_sid.at(i) +=
        result.matching_policy.pairing_offsets_by_3p_sid.at(i - 1);
  }
  result.matching_policy.pairings_by_3p_sid.resize(pairings.size());
  auto next_pairing_offset = result.matching_policy.pairing_offsets_by_3p_sid;
  for (const auto& pairing : pairings) {
    auto& pairing_offset = next_pairing_offset.at(pairing.sid_3p_id);
    result.matching_policy.pairings_by_3p_sid.at(pairing_offset) = pairing.pairing;
    ++pairing_offset;
  }

  result.matching_policy.sample_by_5p_sid.resize(result.sid_5p_pool.size());
  result.matching_policy.sample_by_3p_sid.resize(result.sid_3p_pool.size());
  for (std::size_t sample_id = 0; sample_id < entries.size(); ++sample_id) {
    const auto& entry = entries.at(sample_id);
    if (sid_5p_counts.at(entry.sid_5p) == 1) {
      result.matching_policy.sample_by_5p_sid.at(sid_5p_id_by_sequence.at(entry.sid_5p)) = static_cast<u32>(sample_id);
    }
    if (sid_3p_counts.at(entry.sid_3p) == 1) {
      result.matching_policy.sample_by_3p_sid.at(sid_3p_id_by_sequence.at(entry.sid_3p)) = static_cast<u32>(sample_id);
    }
  }
  return result;
}

BarcodePool LoadSampleSheet(const fs::path& sample_sheet) {
  return LoadSampleSheetSidPools(sample_sheet, DualSidMode::kNone).sid_5p_pool;
}

}  // namespace xoos::demux
