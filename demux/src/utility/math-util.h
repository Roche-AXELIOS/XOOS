#pragma once

namespace xoos::demux {

template <class T>
inline T CeilDiv(T x, T y) {
  return (x + y - 1) / y;
}

template <class T, class U>
inline T RoundUp(T x, U y) {
  return y * ((x + y - 1) / y);
}

}  // namespace xoos::demux
