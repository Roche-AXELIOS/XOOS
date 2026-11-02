#pragma once

#include <cstddef>

namespace xoos::demux::adapter_component_index {

struct Duplex {
  static constexpr std::size_t kRunway = 0U;
  static constexpr std::size_t kSidFwd = 1U;
  static constexpr std::size_t kStart = 2U;
  static constexpr std::size_t kLoop = 3U;
  static constexpr std::size_t kAnchorStop = 4U;
  static constexpr std::size_t kSidRev = 5U;
};

struct DuplexStem {
  static constexpr std::size_t kRunway = 0U;
  static constexpr std::size_t kSidFwd = 1U;
  static constexpr std::size_t kStart = 2U;
  static constexpr std::size_t kStemPre = 3U;
  static constexpr std::size_t kLoop = 4U;
  static constexpr std::size_t kSidRev = 5U;
  static constexpr std::size_t kAnchorStop = 6U;
  static constexpr std::size_t kStemPost = 7U;
};

struct DuplexUmi {
  static constexpr std::size_t kRunway = 0U;
  static constexpr std::size_t kSidFwd = 1U;
  static constexpr std::size_t kStart = 2U;
  static constexpr std::size_t kStemPre = 3U;
  static constexpr std::size_t kLoop = 4U;
  static constexpr std::size_t kUmiFwd = 5U;
  static constexpr std::size_t kSidRev = 6U;
  static constexpr std::size_t kAnchorStop = 7U;
  static constexpr std::size_t kStemPost = 8U;
  static constexpr std::size_t kUmiRev = 9U;
};

struct Ysu {
  static constexpr std::size_t kSidFwd = 0U;
  static constexpr std::size_t kUmiFwd = 1U;
  static constexpr std::size_t kUmiRev = 2U;
  static constexpr std::size_t kSidRev = 3U;
};

struct Ys {
  static constexpr std::size_t kRunway = 0U;
  static constexpr std::size_t kSidFwd = 1U;
  static constexpr std::size_t kSpacer = 2U;
  static constexpr std::size_t kSidRev = 3U;
  static constexpr std::size_t kBait = 4U;
};

struct Simplex {
  static constexpr std::size_t kAnchorFwd = 0U;
  static constexpr std::size_t kSidFwd = 1U;
  static constexpr std::size_t kStemFwd = 2U;
  static constexpr std::size_t kStemRev = 3U;
  static constexpr std::size_t kSidRev = 4U;
  static constexpr std::size_t kAnchorRev = 5U;
};

struct Simplex10x {
  static constexpr std::size_t kAnchor = 0U;
  static constexpr std::size_t kSid = 1U;
  static constexpr std::size_t kStem = 2U;
};

}  // namespace xoos::demux::adapter_component_index
