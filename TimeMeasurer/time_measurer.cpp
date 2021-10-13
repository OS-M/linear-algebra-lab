#include "time_measurer.h"

TimeMeasurer::TimeMeasurer() {
  start = std::chrono::high_resolution_clock::now();
}

std::string TimeMeasurer::GetDurationString() const {
  auto now = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = now - start;
  return "Duration was " + std::to_string(duration.count()) + "s";
}
