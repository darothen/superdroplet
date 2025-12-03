#include "utils/time.hpp"
#include <sstream>
#include <iomanip>

namespace sd_cpp {
namespace utils {

Stopwatch::Stopwatch(unsigned int seconds) {
    const unsigned int seconds_over = seconds % 60;
    minutes = (seconds - seconds_over) / 60;
    this->seconds = seconds_over;
}

unsigned int Stopwatch::total_seconds() const {
    return minutes * 60 + seconds;
}

void Stopwatch::increment(unsigned int seconds_to_add) {
    seconds += seconds_to_add;
    if (seconds >= 60) {
        const unsigned int seconds_overflow = seconds % 60;
        const unsigned int minutes_overflow = (seconds - seconds_overflow) / 60;
        minutes += minutes_overflow;
        seconds = seconds_overflow;
    }
}

std::string Stopwatch::to_string() const {
    std::ostringstream oss;
    oss << std::setw(4) << minutes << " minutes, " 
        << std::setw(2) << seconds << " seconds";
    return oss.str();
}

} // namespace utils
} // namespace sd_cpp
