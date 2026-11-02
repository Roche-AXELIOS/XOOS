#include "adapter-dsl/dsl/dsl-compiler.h"

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <xoos/error/error.h>

#include "adapter-dsl/model/adapter-model-constants.h"
#include "adapter-dsl/model/adapter-model.h"

namespace xoos::demux::dsl {

namespace {

// "sid" / "umi" / "ins"
constexpr std::size_t kDslTypePrefixLength = 3;

/// Maximum number of duplex arms allowed in a DSL string (two `(...)` groups).
constexpr std::size_t kMaxDuplexArms = 2;

/// A minimal cursor over the DSL string with position tracking for error messages.
class Cursor {
 public:
  explicit Cursor(const std::string_view text) : _text(text) {
  }

  bool AtEnd() const {
    return _pos >= _text.size();
  }

  char Peek() const {
    return _text[_pos];
  }

  char Next() {
    const char c = _text[_pos];
    ++_pos;
    return c;
  }

  std::size_t Pos() const {
    return _pos;
  }

  std::string_view Text() const {
    return _text;
  }

 private:
  std::string_view _text;
  std::size_t _pos{0};
};

/**
 * Order-preserving id -> declared entry type table. sid/umi is known only from the DSL bracket,
 * not from the read_structure id token, so we capture it while parsing brackets and use it to
 * build the adapters map (and to enforce that an id is not declared as both sid and umi).
 */
class IdTypeTable {
 public:
  void Record(const std::string& id, const EntryType type) {
    for (const auto& [existing_id, existing_type] : _entries) {
      if (existing_id == id) {
        if (existing_type != type) {
          throw error::Error(
              "DSL id '{}' is declared as both '{}' and '{}'", id, ToString(existing_type), ToString(type));
        }
        return;
      }
    }
    _entries.emplace_back(id, type);
  }

  std::vector<std::pair<std::string, AdapterEntry>> ToBareAdapters() const {
    std::vector<std::pair<std::string, AdapterEntry>> adapters;
    adapters.reserve(_entries.size());
    for (const auto& [id, type] : _entries) {
      AdapterEntry entry;
      entry.type = type;
      adapters.emplace_back(id, std::move(entry));
    }
    return adapters;
  }

 private:
  std::vector<std::pair<std::string, EntryType>> _entries;
};

/**
 * Parse the contents of a `[...]` bracket into a leaf token. The opening '[' has been consumed;
 * this consumes through the matching ']'. sid/umi ids are recorded into @p id_types.
 */
LeafToken ParseBracket(Cursor& cursor, IdTypeTable& id_types) {
  std::string body;
  bool closed = false;
  while (!cursor.AtEnd()) {
    const char c = cursor.Next();
    if (c == kBracketClose) {
      closed = true;
      break;
    }
    body += c;
  }
  if (!closed) {
    throw error::Error("DSL bracket is not closed with ']': '{}'", cursor.Text());
  }
  if (body.size() < kDslTypePrefixLength) {
    throw error::Error("DSL bracket '[{}]' is too short; expected a 3-character type prefix (sid/umi/ins)", body);
  }
  const std::string type_prefix = body.substr(0, kDslTypePrefixLength);
  const std::string rest = body.substr(kDslTypePrefixLength);

  if (type_prefix == kDslInsertKeyword) {
    // Reserved keyword: `[ins]`, `[ins+]`, or `[ins-]`. No id, no colon.
    if (rest.empty()) {
      return InsertToken{std::nullopt};
    }
    if (rest == std::string{kOrientationForward}) {
      return InsertToken{Orientation::kForward};
    }
    if (rest == std::string{kOrientationReverse}) {
      return InsertToken{Orientation::kReverse};
    }
    throw error::Error("DSL insert '[{}]' may only be '[ins]', '[ins+]', or '[ins-]'", body);
  }

  const EntryType entry_type = [&]() {
    if (type_prefix == kDslSid) {
      return EntryType::kSid;
    }
    if (type_prefix == kDslUmi) {
      return EntryType::kUmi;
    }
    throw error::Error(
        "DSL bracket '[{}]' has unknown type prefix '{}' (expected sid, umi, or ins)", body, type_prefix);
  }();

  // sid / umi require `:id` with an optional trailing +/-.
  if (rest.empty() || rest.front() != kColon) {
    throw error::Error("DSL bracket '[{}]' of type '{}' must be followed by ':<id>'", body, type_prefix);
  }
  std::string id = rest.substr(1);
  if (id.empty()) {
    throw error::Error("DSL bracket '[{}]' has an empty id", body);
  }
  Orientation orientation = Orientation::kForward;
  const char last = id.back();
  if (last == kOrientationForward) {
    id.pop_back();
    orientation = Orientation::kForward;
  } else if (last == kOrientationReverse) {
    id.pop_back();
    orientation = Orientation::kReverse;
  }
  if (!IsValidIdBase(id)) {
    throw error::Error("DSL id '{}' in '[{}]' must match ^[a-z][a-z0-9_]+$", id, body);
  }
  id_types.Record(id, entry_type);
  return IdToken{std::move(id), orientation};
}

/// Parse a run of uppercase ACGT characters into a fixed token. At least one has been peeked.
FixedToken ParseFixed(Cursor& cursor) {
  std::string sequence;
  while (!cursor.AtEnd()) {
    const char c = cursor.Peek();
    if (c >= 'A' && c <= 'Z') {
      sequence += cursor.Next();
    } else {
      break;
    }
  }
  if (!IsUpperAcgt(sequence)) {
    throw error::Error("DSL fixed sequence '{}' must be uppercase ACGT", sequence);
  }
  return FixedToken{std::move(sequence)};
}

/// The staged parse output: a flat list of top-level items where `(` / `)` are markers.
struct ParenOpen {};

struct ParenClose {};

using StagedItem = std::variant<LeafToken, ParenOpen, ParenClose>;

/// Lex the DSL into a flat staged item stream (leaves plus paren markers), recording id types.
std::vector<StagedItem> Lex(std::string_view dsl_text, IdTypeTable& id_types) {
  Cursor cursor{dsl_text};
  std::vector<StagedItem> items;
  while (!cursor.AtEnd()) {
    const char c = cursor.Peek();
    if (c == kBracketOpen) {
      cursor.Next();
      items.emplace_back(ParseBracket(cursor, id_types));
    } else if (c == kParenOpen) {
      cursor.Next();
      items.emplace_back(ParenOpen{});
    } else if (c == kParenClose) {
      cursor.Next();
      items.emplace_back(ParenClose{});
    } else if (c >= 'A' && c <= 'Z') {
      items.emplace_back(LeafToken{ParseFixed(cursor)});
    } else {
      throw error::Error("DSL has an unexpected character '{}' at position {} in '{}'", c, cursor.Pos(), dsl_text);
    }
  }
  return items;
}

/// Lift a leaf token into a top-level read-structure Token.
Token LeafToTopToken(const LeafToken& leaf) {
  if (const auto* const fixed = std::get_if<FixedToken>(&leaf)) {
    return Token{*fixed};
  }
  if (const auto* const id = std::get_if<IdToken>(&leaf)) {
    return Token{*id};
  }
  return Token{std::get<InsertToken>(leaf)};
}

/// Assemble the read_structure from the staged item stream, folding `(...)` groups into arms.
std::vector<Token> Assemble(const std::vector<StagedItem>& items, const std::string_view dsl_text) {
  std::vector<Token> read_structure;
  bool in_arm = false;
  std::size_t arm_count = 0;
  ArmToken current_arm;

  for (const StagedItem& item : items) {
    if (std::holds_alternative<ParenOpen>(item)) {
      if (in_arm) {
        throw error::Error("DSL has a nested '(' which is not supported: '{}'", dsl_text);
      }
      ++arm_count;
      // Reachable with 3+ '(' groups (e.g. "(a)(b)(c)"); Sonar's engine misjudges the loop bound.
      if (arm_count > kMaxDuplexArms) {  // NOSONAR(M23_002)
        throw error::Error("DSL has more than two duplex arms (only two '(...)' groups are allowed): '{}'", dsl_text);
      }
      in_arm = true;
      current_arm = ArmToken{};
    } else if (std::holds_alternative<ParenClose>(item)) {
      if (!in_arm) {
        throw error::Error("DSL has an unmatched ')': '{}'", dsl_text);
      }
      if (current_arm.tokens.empty()) {
        throw error::Error("DSL duplex arm '()' is empty: '{}'", dsl_text);
      }
      read_structure.emplace_back(std::move(current_arm));
      current_arm = ArmToken{};
      in_arm = false;
    } else {
      const auto& leaf = std::get<LeafToken>(item);
      if (in_arm) {
        current_arm.tokens.push_back(leaf);
      } else {
        read_structure.emplace_back(LeafToTopToken(leaf));
      }
    }
  }
  if (in_arm) {
    throw error::Error("DSL has an unmatched '(' (missing closing ')'): '{}'", dsl_text);
  }
  if (arm_count == 1) {
    throw error::Error("DSL has a single duplex arm; a duplex fold needs two '(...)' arms: '{}'", dsl_text);
  }
  return read_structure;
}

}  // namespace

CompiledAdapter CompileDsl(const std::string_view dsl_text) {
  if (dsl_text.empty()) {
    throw error::Error("DSL string is empty");
  }
  IdTypeTable id_types;
  const std::vector<StagedItem> items = Lex(dsl_text, id_types);
  CompiledAdapter adapter;
  adapter.read_structure = Assemble(items, dsl_text);
  adapter.adapters = id_types.ToBareAdapters();
  return adapter;
}

AdapterBundle CompileDslToBundle(const std::string& name, const std::string_view dsl_text) {
  AdapterBundle bundle;
  bundle.adapters.emplace_back(name, CompileDsl(dsl_text));
  ValidateBundle(bundle);
  return bundle;
}

}  // namespace xoos::demux::dsl
