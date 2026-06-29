#pragma once

#include <chrono>

namespace xoos::demux {

/**
 * A simple stopwatch class for measuring elapsed time.
 */
class StopWatch {
 private:
  using Clock = std::chrono::high_resolution_clock;

 public:
  StopWatch() : _start_point{Clock::now()} {}

  /**
   * @brief Calculates and returns the elapsed time in the specified duration unit.
   *
   * This function calculates the duration between the current time and the start point
   * of the stopwatch. The calculated duration is then converted to the specified duration
   * unit using duration_cast and the count (numeric value) of the elapsed time is returned.
   *
   * @tparam U The duration unit type. Defaults to nanoseconds.
   * @return The elapsed time in the specified duration unit.
   */
  template <class U = std::chrono::nanoseconds>
  Clock::duration::rep ElapsedTime() const {
    return std::chrono::duration_cast<U>(Clock::now() - _start_point).count();
  }

 private:
  Clock::time_point _start_point;
};

}  // namespace xoos::demux
