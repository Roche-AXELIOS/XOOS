#pragma once

#include <string>

#include "adapter-design/adapter-design.h"
#include "adapter-dsl/model/adapter-model.h"

namespace xoos::demux::dsl {

/**
 * @brief Assemble a runtime @ref AdapterDesign blueprint from a compiled DSL adapter.
 *
 * This is the config loader's final hop: it turns a validated @ref CompiledAdapter into the flat
 * @ref AdapterDesign the existing bundle loader consumes. Like the manifest it replaces, the result
 * is a pure blueprint -- @ref AdapterComponentDefinition entries carry FASTA *paths* or fixed
 * literals, not loaded pools. The bundle loader (@ref LoadAdapterDesignBundle) still loads every
 * pool once, unchanged, and still applies the sample-sheet SID override.
 *
 * Variable matcher components are classified as SID or UMI. Fixed literals are a separate component
 * type and are consumed by position, not by a sub-role tag.
 *
 * Geometry (which components land on the 5'/3' side, the duplex anchor/stem/loop split) is read
 * arm-aware from @ref CompiledAdapter::read_structure. Derived values -- UMI flank affixes and YSU
 * sid baits -- are taken from a single derivation site: @ref ResolveStructure followed by
 * @ref Infer. Nothing is derived a second time here.
 *
 * @throws xoos::Error on a structural violation (strict duplex arm symmetry is enforced here, at
 * selection time) or an unhandled design type.
 */
AdapterDesign ToAdapterDesign(const std::string& name, const CompiledAdapter& adapter);

}  // namespace xoos::demux::dsl
