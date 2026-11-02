#pragma once

#include <string>
#include <string_view>
#include <utility>

namespace xoos::demux {

/**
 * @brief A single fixed adapter sequence (anchor, stem, runway, loop).
 *
 * Fixed components are constant literals resolved by the DSL assembler. Unlike SeqMatcher, no hash
 * table or edit-distance expansion is built -- the sequence is stored as a plain string and
 * consumed directly by the trimmers. The literal is already in molecule read order, so no
 * orientation, transform, or search metadata is retained here.
 */
class FixedSequence {
 public:
  FixedSequence() = default;

  explicit FixedSequence(std::string sequence) : _sequence{std::move(sequence)} {
  }

  std::string_view Seq() const {
    return _sequence;
  }

  size_t Size() const {
    return _sequence.size();
  }

 private:
  std::string _sequence;
};

}  // namespace xoos::demux
