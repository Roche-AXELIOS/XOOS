#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <utility>

#include <xoos/types/vec.h>
#include <xoos/vfs/vfs.h>

#include "adapter-design.h"
#include "sequence/fixed-sequence.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

using VfsPtr = vfs::VirtualFilesystemPtr;

/**
 * @brief One adapter component in molecule read order, in the form the trimmers consume.
 *
 * A component is either a variable barcode (sid / umi), which carries a @ref matcher (a LUT), or a
 * fixed literal, which carries a @ref fixed sequence. The two are mutually exclusive and the
 * distinction is the matcher: @ref matcher is non-null for a variable component and null for a fixed
 * one -- no separate type tag is stored because no consumer branches on it (every trimmer addresses
 * a known position and asks for exactly the form it expects). @ref bait holds a sid's Infer'd bait
 * seed when present. The trimmers address components by their position in the vector; there is no
 * 5'/3' side tag, and orientation is already baked into each component by the assembler.
 */
struct AdapterComponent {
  /// @brief The fixed literal, valid iff @ref matcher is null.
  FixedSequence fixed;
  /// @brief The LUT matcher, non-null iff this is a variable (sid / umi) component.
  SeqMatcherPtr matcher;
  /// @brief The Infer'd bait seed for a sid component; empty when the component has none.
  std::string bait;

  /// @brief Whether this is a variable component (has a LUT matcher) rather than a fixed literal.
  bool HasMatcher() const {  // NOSONAR(S5536) -- exercised by tests (not in sonar.sources)
    return matcher != nullptr;
  }
};

/**
 * @brief The fixed literal of the component at read-order position @p index within @p components.
 * Throws if the position is out of range or the component carries a matcher (i.e. is not a fixed
 * literal). Free function so trimmers can address a plain component vector without an owning bundle.
 */
const FixedSequence& FixedAt(const vec<AdapterComponent>& components, std::size_t index);

/**
 * @brief The matcher of the component at read-order position @p index within @p components.
 * Throws if the position is out of range or the component has no matcher (i.e. is a fixed literal).
 * Free function so trimmers can address a plain component vector without an owning bundle.
 */
SeqMatcherPtr MatcherAt(const vec<AdapterComponent>& components, std::size_t index);

/**
 * @brief The barcode pool of the first variable (matcher-bearing) component in read order.
 *
 * Used when a caller needs a design-agnostic matcher pool without depending on any one adapter's
 * component layout. Throws if @p components has no variable component.
 */
const BarcodePool& FirstMatcherPool(const vec<AdapterComponent>& components);

/**
 * @brief Try to load the untransformed SID pool referenced by @p design.
 *
 * The returned sequences preserve their canonical FASTA orientation. Matcher-specific transforms
 * are not applied. Returns std::nullopt when the referenced FASTA is unavailable, which allows a
 * supplied sample-sheet pool to override a missing bundle default.
 */
std::optional<BarcodePool> TryLoadCanonicalSidPool(const fs::path& bundle, const AdapterDesign& design);

/**
 * Contains the components (matchers and fixed literals) required to demux and trim a given
 * adapter design, in a single read-order vector.
 */
struct AdapterDesignBundle {
  /// @brief All components in molecule read order (@ref AdapterComponent). Consumed by position.
  vec<AdapterComponent> components;
  AdapterDesign design;

  /**
   * @brief The matcher of the component at read-order position @p index.
   * Throws if the position is out of range or the component carries no matcher (a fixed literal).
   */
  SeqMatcherPtr MatcherAt(std::size_t index) const;  // NOSONAR(S5536) -- used by adapter trimmers

  /**
   * @brief The fixed literal of the component at read-order position @p index.
   * Throws if the position is out of range or the component is not a fixed literal.
   */
  const FixedSequence& FixedAt(std::size_t index) const;
};

/// @brief Overrides one SID component's FASTA-backed pool with an in-memory sample-sheet pool.
struct ComponentBarcodePoolOverride {
  ComponentBarcodePoolOverride(const std::size_t index, BarcodePool barcode_pool)
      : component_index{index}, pool{std::move(barcode_pool)} {
  }

  std::size_t component_index{0};
  BarcodePool pool;
};

using ComponentBarcodePoolOverrides = vec<ComponentBarcodePoolOverride>;

/**
 * Load the @p design from the provided @p bundle, this involves reading in the sequences referenced by the @p design
 * and generating the appropriate matchers according to the definitions contained in the @p design. If an @p sid_pool is
 * provided, all SID matchers will use @p sid_pool instead of the FASTA paths defined in the @p design.
 * SID matcher definitions with empty FASTA paths must be supplied through an override.
 */
AdapterDesignBundle LoadAdapterDesignBundle(const fs::path& bundle,
                                            const AdapterDesign& design,
                                            const std::optional<BarcodePool>& sid_pool,
                                            size_t threads);

/**
 * Load the @p design using per-component SID pool overrides. Overrides are keyed by read-order component index.
 * This is the path used for dual-SID designs whose two sides are supplied by the sample sheet.
 */
AdapterDesignBundle LoadAdapterDesignBundle(const fs::path& bundle,
                                            const AdapterDesign& design,
                                            const ComponentBarcodePoolOverrides& sid_pool_overrides,
                                            size_t threads);

/**
 * Similar to LoadAdapterDesignBundle(const fs::path& bundle, const AdapterDesign& design, size_t threads)
 * but loads the bundle from a virtual filesystem.
 */
AdapterDesignBundle LoadAdapterDesignBundle(const VfsPtr& vfs,
                                            const AdapterDesign& design,
                                            const std::optional<BarcodePool>& sid_pool,
                                            size_t threads);

AdapterDesignBundle LoadAdapterDesignBundle(const VfsPtr& vfs,
                                            const AdapterDesign& design,
                                            const ComponentBarcodePoolOverrides& sid_pool_overrides,
                                            size_t threads);

}  // namespace xoos::demux
