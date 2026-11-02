#pragma once

#include <string>
#include <string_view>

#include <xoos/types/int.h>
#include <xoos/types/str-container.h>
#include <xoos/types/vec.h>

namespace xoos::alignment_metrics {

/**
 * Maximum number of distinct region labels supported. Label membership is stored as a `u16` bit
 * vector (one bit per label), so at most 16 labels can be tracked.
 */
constexpr size_t kMaxLabels = 16;

/**
 * Maps region-label names to bit positions in a `u16` bit vector and back.
 *
 * Labels are assigned bit positions in first-seen order as they are added while parsing a
 * labeled BED file. The resulting bit vector for a region is the union of the bits of every label
 * attached to that region, allowing a genomic position to belong to multiple label classes at once.
 *
 * An empty set (no labels added) represents the non-stratified mode, in which metrics are
 * written to the top-level output directory without any per-label sub-hierarchy.
 */
class RegionLabelsSet {
 public:
  /**
   * Add a label name, returning the single-bit mask (`1u << index`) for its assigned bit
   * position. Adding an already-known label returns its existing mask without allocating a new
   * bit.
   *
   * @param name Label name; must be non-empty and safe to use as a directory name.
   * @return The single-bit `u16` mask for the label.
   * @throws error::Error if the name contains characters that are illegal in a directory name, or if
   *         registering it would exceed `kMaxLabels` distinct labels.
   */
  u16 AddLabel(std::string_view name);

  /**
   * Get the single-bit mask for an already-registered label, or 0 if the label is unknown.
   */
  u16 GetMask(std::string_view name) const;

  /**
   * Number of distinct labels registered.
   */
  size_t Size() const;

  /**
   * True when no labels have been registered (non-stratified mode).
   */
  bool Empty() const;

  /**
   * Get the label name assigned to the given bit index (0-based, in registration order).
   * @throws error::Error if `index` is out of range.
   */
  const std::string& NameAt(size_t index) const;

  /**
   * The label names in registration (bit-index) order.
   */
  const vec<std::string>& Names() const;

  /**
   * Validate that a label name is usable as a directory name. Rejects empty names, path separators,
   * NUL and control characters, any whitespace, and the reserved names "." and "..".
   * @throws error::Error naming the offending label.
   */
  static void ValidateLabelName(std::string_view name);

 private:
  // Label names indexed by assigned bit position (registration order).
  vec<std::string> _names;
  // Label name -> single-bit mask lookup.
  StrUnorderedMap<u16> _mask_by_name;
};

}  // namespace xoos::alignment_metrics
