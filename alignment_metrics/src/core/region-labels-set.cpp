#include "core/region-labels-set.h"

#include <locale>
#include <string>
#include <string_view>

#include <xoos/error/error.h>
#include <xoos/types/int.h>

namespace xoos::alignment_metrics {

void RegionLabelsSet::ValidateLabelName(const std::string_view name) {
  if (name.empty()) {
    throw error::Error("Region label must not be empty");
  }
  if (name == "." || name == "..") {
    throw error::Error("Region label '{}' is reserved and cannot be used as a directory name", name);
  }
  for (const char character : name) {
    if (character == '/' || character == '\\') {
      throw error::Error("Region label '{}' must not contain a path separator ('/' or '\\')", name);
    }
    if (character == '\0' || std::iscntrl(character, std::locale::classic())) {
      throw error::Error("Region label '{}' must not contain control characters", name);
    }
    if (std::isspace(character, std::locale::classic())) {
      throw error::Error("Region label '{}' must not contain whitespace", name);
    }
  }
}

u16 RegionLabelsSet::AddLabel(const std::string_view name) {
  const auto existing = _mask_by_name.find(name);
  if (existing != _mask_by_name.end()) {
    return existing->second;
  }
  ValidateLabelName(name);
  if (_names.size() >= kMaxLabels) {
    throw error::Error(
        "Too many distinct region labels: at most {} are supported, but an additional label '{}' was found",
        kMaxLabels,
        name);
  }
  const auto index = _names.size();
  const auto mask = static_cast<u16>(1U << index);
  _names.emplace_back(name);
  _mask_by_name.try_emplace(std::string(name), mask);
  return mask;
}

u16 RegionLabelsSet::GetMask(const std::string_view name) const {
  const auto existing = _mask_by_name.find(name);
  return existing != _mask_by_name.end() ? existing->second : static_cast<u16>(0);
}

size_t RegionLabelsSet::Size() const {
  return _names.size();
}

bool RegionLabelsSet::Empty() const {
  return _names.empty();
}

const std::string& RegionLabelsSet::NameAt(const size_t index) const {
  if (index >= _names.size()) {
    throw error::Error("Label index {} is out of range ({} labels registered)", index, _names.size());
  }
  return _names.at(index);
}

const vec<std::string>& RegionLabelsSet::Names() const {
  return _names;
}

}  // namespace xoos::alignment_metrics
