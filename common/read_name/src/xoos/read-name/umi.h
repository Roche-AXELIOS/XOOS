#pragma once

#include <cstddef>
#include <optional>
#include <string>

namespace xoos::read_name {

/// A UMI (Unique Molecular Identifier) sequence, or `std::nullopt` when absent.
using Umi = std::optional<std::string>;

/// Represents no UMI in the read name.
const std::string kNoUmi = "*";

struct UmiPair {
  Umi umi5p;
  Umi umi3p;

  auto operator<=>(const UmiPair&) const = default;
};

struct UmiPairHash {
  size_t operator()(const UmiPair& umi_pair) const;
};

/**
 * Check if the UMI is a sequence containing only ACGT characters.
 * @param umi The UMI to check.
 */
bool IsSequenceUmi(const std::optional<std::string>& umi);

}  // namespace xoos::read_name
