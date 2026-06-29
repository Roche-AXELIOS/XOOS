#pragma once

#include <xoos/vfs/vfs.h>

#include <tuple>
#include <vector>

#include "adapter-design.h"
#include "sequence/matcher/seq-matcher.h"

namespace xoos::demux {

using VfsPtr = vfs::VirtualFilesystemPtr;

/**
 * Contains the LUTs required to demux and trim a given
 * adapter design.
 */
struct AdapterDesignBundle {
  std::vector<std::tuple<BarcodeType, SeqMatcherPtr>> adapter_5p;
  SeqMatcherPtr GetMatcher5p(const BarcodeType& type) const;

  std::vector<std::tuple<BarcodeType, SeqMatcherPtr>> adapter_3p;
  SeqMatcherPtr GetMatcher3p(const BarcodeType& type) const;
};

/**
 * Load the @p design from the provided @p bundle, this involves reading in the sequences referenced by the @p design
 * and generating the appropriate LUT according to the definitions contained in the @p design. If an @p sid_pool is
 * provided, than all kSID barcodes will use @p sid_pool instead of those defined in the @p design.
 */
AdapterDesignBundle LoadAdapterDesignBundle(const fs::path& bundle, const AdapterDesign& design,
                                            const std::optional<BarcodePool>& sid_pool, size_t threads);

/**
 * Similar to LoadAdapterDesignBundle(const fs::path& bundle, const AdapterDesign& design, size_t threads)
 * but loads the bundle from a virtual filesystem.
 */
AdapterDesignBundle LoadAdapterDesignBundle(const VfsPtr& vfs, const AdapterDesign& design,
                                            const std::optional<BarcodePool>& sid_pool, size_t threads);

}  // namespace xoos::demux
