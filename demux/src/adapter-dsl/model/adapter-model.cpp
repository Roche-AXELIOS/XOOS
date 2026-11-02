#include "adapter-dsl/model/adapter-model.h"

#include <algorithm>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <xoos/error/error.h>
#include <xoos/types/int.h>
#include <xoos/util/sequence-functions.h>

#include "adapter-dsl/model/adapter-model-constants.h"

namespace xoos::demux::dsl {

// --- Enum string helpers ---

std::string_view ToString(const EntryType type) {
  using enum EntryType;
  switch (type) {
    case kSid:
      return kEntryTypeSid;
    case kUmi:
      return kEntryTypeUmi;
    case kFixed:
      return kEntryTypeFixed;
    default:
      throw error::Error("Unhandled EntryType value {}", std::to_underlying(type));
  }
}

EntryType ParseEntryType(std::string_view str) {
  using enum EntryType;
  if (str == kEntryTypeSid) {
    return kSid;
  }
  if (str == kEntryTypeUmi) {
    return kUmi;
  }
  if (str == kEntryTypeFixed) {
    return kFixed;
  }
  throw error::Error("Unknown adapters entry type '{}', valid values: 'sid', 'umi', 'fixed'", str);
}

std::string_view ToString(const DesignType type) {
  using enum DesignType;
  switch (type) {
    case kDuplex:
      return kDesignTypeDuplex;
    case kDuplexStem:
      return kDesignTypeDuplexStem;
    case kDuplexUmi:
      return kDesignTypeDuplexUmi;
    case kYsu:
      return kDesignTypeYsu;
    case kYs:
      return kDesignTypeYs;
    case kSimplex:
      return kDesignTypeSimplex;
    case kSimplex10x:
      return kDesignTypeSimplex10x;
    default:
      throw error::Error("Unhandled DesignType value {}", std::to_underlying(type));
  }
}

DesignType ParseDesignType(const std::string_view str) {
  using enum DesignType;
  if (str == kDesignTypeDuplex) {
    return kDuplex;
  }
  if (str == kDesignTypeDuplexStem) {
    return kDuplexStem;
  }
  if (str == kDesignTypeDuplexUmi) {
    return kDuplexUmi;
  }
  if (str == kDesignTypeYsu) {
    return kYsu;
  }
  if (str == kDesignTypeYs) {
    return kYs;
  }
  if (str == kDesignTypeSimplex) {
    return kSimplex;
  }
  if (str == kDesignTypeSimplex10x) {
    return kSimplex10x;
  }
  throw error::Error(
      "Unknown adapter design type '{}', valid values: "
      "'{}', '{}', '{}', '{}', '{}', '{}', '{}'",
      str,
      kDesignTypeDuplex,
      kDesignTypeDuplexStem,
      kDesignTypeDuplexUmi,
      kDesignTypeYsu,
      kDesignTypeYs,
      kDesignTypeSimplex,
      kDesignTypeSimplex10x);
}

char OrientationSuffix(const Orientation orientation) {
  return orientation == Orientation::kForward ? kOrientationForward : kOrientationReverse;
}

bool IsUpperAcgt(const std::string_view s) {
  return !s.empty() &&
         std::ranges::all_of(s, [](const char c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T'; });
}

/// Minimum length of an id base: one leading letter plus at least one more character.
constexpr std::size_t kMinIdBaseLength = 2;

namespace {

bool IsIdBaseChar(const char c) {
  const bool is_lower = c >= 'a' && c <= 'z';
  const bool is_digit = c >= '0' && c <= '9';
  return is_lower || is_digit || c == '_';
}

}  // namespace

bool IsValidIdBase(const std::string_view s) {
  // ^[a-z][a-z0-9_]+$ : at least two characters, first is a-z, rest a-z/0-9/_.
  if (s.size() < kMinIdBaseLength) {
    return false;
  }
  if (s.front() < 'a' || s.front() > 'z') {
    return false;
  }
  return std::ranges::all_of(s.substr(1), IsIdBaseChar);
}

// --- Leaf token <-> compiled string ---

std::string LeafToCompiledString(const LeafToken& leaf) {
  if (const auto* const fixed = std::get_if<FixedToken>(&leaf)) {
    return fixed->sequence;
  }
  if (const auto* const id = std::get_if<IdToken>(&leaf)) {
    return id->base + OrientationSuffix(id->orientation);
  }
  const auto& ins = std::get<InsertToken>(leaf);
  std::string out{kInsertKeyword};
  if (ins.orientation.has_value()) {
    out += OrientationSuffix(ins.orientation.value());
  }
  return out;
}

namespace {

/**
 * Split a trailing +/- orientation marker off a compiled string token.
 * Returns the base (without the marker) and the parsed orientation, if any.
 */
std::pair<std::string_view, std::optional<Orientation>> SplitOrientation(std::string_view token) {
  if (token.empty()) {
    return {token, std::nullopt};
  }
  const char last = token.back();
  if (last == kOrientationForward) {
    return {token.substr(0, token.size() - 1), Orientation::kForward};
  }
  if (last == kOrientationReverse) {
    return {token.substr(0, token.size() - 1), Orientation::kReverse};
  }
  return {token, std::nullopt};
}

}  // namespace

LeafToken LeafFromCompiledString(std::string_view token) {
  if (token.empty()) {
    throw error::Error("Empty read_structure token");
  }
  // Case = class: an uppercase ACGT run is a fixed literal (the data itself).
  if (IsUpperAcgt(token)) {
    return FixedToken{std::string{token}};
  }
  const auto [base, orientation] = SplitOrientation(token);
  if (base == kInsertKeyword) {
    return InsertToken{orientation};
  }
  if (!IsValidIdBase(base)) {
    throw error::Error(
        "read_structure token '{}' is neither uppercase ACGT, the insert keyword, nor a valid id "
        "(^[a-z][a-z0-9_]+$)",
        token);
  }
  // A bare id with no marker is forward (+ and no marker are equivalent).
  return IdToken{std::string{base}, orientation.value_or(Orientation::kForward)};
}

// --- CompiledAdapter helpers ---

const AdapterEntry* CompiledAdapter::Find(const std::string_view base) const {
  for (const auto& [name, entry] : adapters) {
    if (name == base) {
      return &entry;
    }
  }
  return nullptr;
}

bool CompiledAdapter::IsDuplex() const {
  return std::ranges::any_of(read_structure, [](const Token& t) { return std::holds_alternative<ArmToken>(t); });
}

// --- AdapterBundle lookup ---

const CompiledAdapter& AdapterBundle::At(const std::string& name) const {
  for (const auto& [adapter_name, adapter] : adapters) {
    if (adapter_name == name) {
      return adapter;
    }
  }
  throw error::Error("Adapter '{}' does not exist in bundle", name);
}

bool AdapterBundle::Contains(const std::string& name) const {
  return std::ranges::any_of(adapters, [&name](const auto& kv) { return kv.first == name; });
}

// --- Validation ---

namespace {

/**
 * A duplex read structure has exactly two arms separated by a single fixed loop literal, so the
 * arm token positions differ by two (arm, loop, arm).
 */
constexpr std::size_t kDuplexArmCount = 2;
constexpr std::size_t kArmLoopArmGap = 2;

/// Visit every leaf token in a read structure in order (descending into arms), calling `fn`.
template <typename Fn>
void ForEachLeaf(const std::vector<Token>& read_structure, const Fn& fn) {
  for (const Token& token : read_structure) {
    std::visit(
        [&fn]<typename T>(const T& value) {
          if constexpr (std::is_same_v<T, ArmToken>) {
            for (const LeafToken& leaf : value.tokens) {
              fn(leaf);
            }
          } else {
            fn(LeafToken{value});
          }
        },
        token);
  }
}

/// Validate a single leaf token: case=class rule, id resolution, id-base regex.
void ValidateLeaf(const std::string& name, const LeafToken& leaf, const CompiledAdapter& adapter) {
  if (const auto* const fixed = std::get_if<FixedToken>(&leaf)) {
    if (!IsUpperAcgt(fixed->sequence)) {
      throw error::Error("Adapter '{}': fixed token '{}' must be uppercase ACGT", name, fixed->sequence);
    }
    return;
  }
  if (const auto* const id = std::get_if<IdToken>(&leaf)) {
    if (!IsValidIdBase(id->base)) {
      throw error::Error("Adapter '{}': id '{}' must match ^[a-z][a-z0-9_]+$", name, id->base);
    }
    if (adapter.Find(id->base) == nullptr) {
      throw error::Error("Adapter '{}': id '{}' has no entry in adapters", name, id->base);
    }
    return;
  }
  // InsertToken: nothing to resolve (self-describing keyword).
}

/// Count how many insert leaf tokens the read structure contains.
u32 CountInserts(const std::vector<LeafToken>& leaves) {
  u32 count = 0;
  for (const LeafToken& leaf : leaves) {
    if (std::holds_alternative<InsertToken>(leaf)) {
      ++count;
    }
  }
  return count;
}

/// Count how many insert leaf tokens the read structure contains.
u32 CountInserts(const CompiledAdapter& adapter) {
  u32 count = 0;
  ForEachLeaf(adapter.read_structure, [&count](const LeafToken& leaf) {
    if (std::holds_alternative<InsertToken>(leaf)) {
      ++count;
    }
  });
  return count;
}

/// Validate the spec's insert-count contract: simplex has one insert; duplex has one insert per arm.
void ValidateInserts(const std::string& name, const CompiledAdapter& adapter) {
  if (!adapter.IsDuplex()) {
    const u32 insert_count = CountInserts(adapter);
    if (insert_count != 1) {
      throw error::Error(
          "Adapter '{}': simplex read_structure must contain exactly one insert, found {}", name, insert_count);
    }
    return;
  }

  for (const Token& token : adapter.read_structure) {
    if (std::holds_alternative<InsertToken>(token)) {
      throw error::Error("Adapter '{}': duplex insert tokens must appear inside duplex arms", name);
    }
  }

  std::size_t arm_index = 0;
  for (const Token& token : adapter.read_structure) {
    if (const auto* const arm = std::get_if<ArmToken>(&token)) {
      const u32 insert_count = CountInserts(arm->tokens);
      if (insert_count != 1) {
        throw error::Error(
            "Adapter '{}': duplex arm {} must contain exactly one insert, found {}", name, arm_index + 1, insert_count);
      }
      ++arm_index;
    }
  }
}

/// Extract the ordered leaf tokens of a single arm.
const std::vector<LeafToken>& ArmTokens(const Token& token) {
  return std::get<ArmToken>(token).tokens;
}

/// A mirrored fixed pair must be reverse complements of each other.
void ValidateMirroredFixed(const std::string& name, const FixedToken& fa, const FixedToken& fb) {
  if (sequence::ReverseComplement(fa.sequence) != fb.sequence) {
    throw error::Error(
        "Adapter '{}': duplex arm fixed sequences '{}' and '{}' are not reverse complements of each other",
        name,
        fa.sequence,
        fb.sequence);
  }
}

/// A mirrored id pair must share a base name and appear on opposite strands.
void ValidateMirroredId(const std::string& name, const IdToken& ia, const IdToken& ib) {
  if (ia.base != ib.base) {
    throw error::Error("Adapter '{}': duplex arm ids '{}' and '{}' do not pair by base name", name, ia.base, ib.base);
  }
  if (ia.orientation == ib.orientation) {
    throw error::Error(
        "Adapter '{}': duplex arm id '{}' must appear on opposite strands in the two arms", name, ia.base);
  }
}

/**
 * Validate a single mirrored pair of leaf tokens from the two duplex arms. `a` and `b` are known
 * to hold the same alternative (their variant indices were compared by the caller), so the
 * `std::get<T>(b)` lookups below are safe. InsertToken pairs positionally with no further constraint.
 */
void ValidateMirroredPair(const std::string& name, const LeafToken& a, const LeafToken& b) {
  std::visit(
      [&name, &b]<typename T>([[maybe_unused]] const T& va) {
        if constexpr (std::is_same_v<T, FixedToken>) {
          ValidateMirroredFixed(name, va, std::get<FixedToken>(b));
        } else if constexpr (std::is_same_v<T, IdToken>) {
          ValidateMirroredId(name, va, std::get<IdToken>(b));
        }
      },
      a);
}

/**
 * The two duplex arms must mirror: reading arm1 left->right and arm2 right->left, each fixed
 * literal in arm1 must be the reverse complement of the corresponding fixed literal in arm2, and
 * each id base in arm1 must equal the id base at the mirrored position in arm2 (opposite strand).
 */
void ValidateArmSymmetry(const std::string& name,
                         const std::vector<LeafToken>& arm1,
                         const std::vector<LeafToken>& arm2) {
  if (arm1.size() != arm2.size()) {
    throw error::Error("Adapter '{}': duplex arms have different lengths ({} vs {}); they must be mirror images",
                       name,
                       arm1.size(),
                       arm2.size());
  }
  const std::size_t n = arm1.size();
  for (std::size_t i = 0; i < n; ++i) {
    const LeafToken& a = arm1[i];
    // mirrored position
    const LeafToken& b = arm2[n - 1 - i];
    if (a.index() != b.index()) {
      throw error::Error("Adapter '{}': duplex arms are not mirror images at paired position {}", name, i);
    }
    ValidateMirroredPair(name, a, b);
  }
}

/**
 * Validate the duplex fold: exactly two arms, mirrored, with a fixed loop literal between them.
 * @param check_arm_symmetry When false, the structural checks (arm count, loop-by-position) still
 * run but the strict reverse-complement arm-symmetry check is skipped. This lets a bundle load a
 * design with asymmetric arms; the strict check is re-applied at selection time.
 */
void ValidateDuplex(const std::string& name, const CompiledAdapter& adapter, const bool check_arm_symmetry) {
  std::vector<std::size_t> arm_positions;
  for (std::size_t i = 0; i < adapter.read_structure.size(); ++i) {
    if (std::holds_alternative<ArmToken>(adapter.read_structure[i])) {
      arm_positions.emplace_back(i);
    }
  }
  if (arm_positions.empty()) {
    return;  // simplex
  }
  if (arm_positions.size() != kDuplexArmCount) {
    throw error::Error("Adapter '{}': duplex must have exactly two arms, found {}", name, arm_positions.size());
  }
  // The loop is the fixed literal between the two arms, recovered by position: arm, loop, arm.
  const std::size_t gap = arm_positions[1] - arm_positions[0];
  if (gap != kArmLoopArmGap) {
    throw error::Error("Adapter '{}': duplex arms must be separated by exactly one fixed loop literal", name);
  }
  const Token& between = adapter.read_structure[arm_positions[0] + 1];
  if (!std::holds_alternative<FixedToken>(between)) {
    throw error::Error("Adapter '{}': the element between duplex arms (the loop) must be a fixed literal", name);
  }
  if (check_arm_symmetry) {
    ValidateArmSymmetry(
        name, ArmTokens(adapter.read_structure[arm_positions[0]]), ArmTokens(adapter.read_structure[arm_positions[1]]));
  }
}

/**
 * True when a config entry carries any v1 matcher knob. Fixed side-cars are metadata-only and must
 * not set these fields.
 */
bool HasMatcherTuning(const AdapterEntry& entry) {
  if (entry.max_edit_distance.has_value()) {
    return true;
  }
  if (entry.wiggle_forward.has_value()) {
    return true;
  }
  if (entry.wiggle_reverse.has_value()) {
    return true;
  }
  if (entry.bait_length_forward.has_value()) {
    return true;
  }
  if (entry.bait_length_reverse.has_value()) {
    return true;
  }
  if (entry.additional_prefix_length.has_value()) {
    return true;
  }
  return entry.additional_suffix_length.has_value();
}

/**
 * Every id base that appears must have an entry, and every entry must be anchored in the read
 * structure: sid/umi entries by id token, fixed side-cars by exact fixed literal.
 */
void ValidateLinkage(const std::string& name, const CompiledAdapter& adapter) {
  std::set<std::string, std::less<>> referenced_ids;
  std::set<std::string, std::less<>> fixed_sequences;
  ForEachLeaf(adapter.read_structure, [&referenced_ids, &fixed_sequences](const LeafToken& leaf) {
    std::visit(
        [&referenced_ids, &fixed_sequences]<typename T>([[maybe_unused]] const T& value) {
          if constexpr (std::is_same_v<T, IdToken>) {
            referenced_ids.insert(value.base);
          } else if constexpr (std::is_same_v<T, FixedToken>) {
            fixed_sequences.insert(value.sequence);
          }
        },
        leaf);
  });
  for (const auto& [key, entry] : adapter.adapters) {
    if (entry.type == EntryType::kFixed) {
      if (!fixed_sequences.contains(key)) {
        throw error::Error(
            "Adapter '{}': fixed adapters entry '{}' does not match any fixed literal in the read structure",
            name,
            key);
      }
      if (entry.sequence != key) {
        throw error::Error(
            "Adapter '{}': fixed adapters entry '{}' sequence '{}' must match its key", name, key, entry.sequence);
      }
      if (HasMatcherTuning(entry)) {
        throw error::Error("Adapter '{}': fixed adapters entry '{}' must not contain matcher tuning fields", name, key);
      }
    } else if (!referenced_ids.contains(key)) {
      throw error::Error("Adapter '{}': adapters entry '{}' is not referenced in the read structure", name, key);
    }
  }
}

}  // namespace

void ValidateCompiledAdapter(const std::string& name, const CompiledAdapter& adapter, const bool check_arm_symmetry) {
  if (adapter.read_structure.empty()) {
    throw error::Error("Adapter '{}': read_structure must not be empty", name);
  }
  ForEachLeaf(adapter.read_structure, [&name, &adapter](const LeafToken& leaf) { ValidateLeaf(name, leaf, adapter); });
  ValidateDuplex(name, adapter, check_arm_symmetry);
  ValidateInserts(name, adapter);
  ValidateLinkage(name, adapter);
}

void ValidateCompiledAdapter(const std::string& name, const CompiledAdapter& adapter) {
  ValidateCompiledAdapter(name, adapter, /*check_arm_symmetry=*/true);
}

void ValidateBundle(const AdapterBundle& bundle) {
  if (bundle.adapters.empty()) {
    throw error::Error("Adapter bundle is empty");
  }
  std::set<std::string, std::less<>> names;
  for (const auto& [name, adapter] : bundle.adapters) {
    if (name.empty()) {
      throw error::Error("Adapter bundle contains an adapter with an empty name");
    }
    if (!names.insert(name).second) {
      throw error::Error("Adapter bundle contains duplicate adapter name '{}'", name);
    }
    // Structural validation only at load time. Strict arm symmetry is deferred to selection time
    // (ToAdapterDesign) so a bundle may legitimately contain designs with asymmetric arms.
    ValidateCompiledAdapter(name, adapter, /*check_arm_symmetry=*/false);
  }
}

}  // namespace xoos::demux::dsl
