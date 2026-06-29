#pragma once

#include "adapter-design/adapter-design-bundle.h"

namespace xoos::demux {

template <typename LutBundle>
LutBundle CreateLutBundle(const AdapterDesignBundle& designs);

template <class LutBundle>
LutBundle LoadLutBundle(const fs::path& bundle, const AdapterDesign& design, const std::optional<BarcodePool>& sid_pool,
                        const size_t threads) {
  const auto adapter_designs = LoadAdapterDesignBundle(bundle, design, sid_pool, threads);
  return CreateLutBundle<LutBundle>(adapter_designs);
}

}  // namespace xoos::demux
