#pragma once

#include <chrono>
#include <string>

class TimeMeasurer {
 public:
  TimeMeasurer();
  std::string GetDurationString() const;

 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start;
};
