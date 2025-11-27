#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <string>

namespace Stats {

// Global verbose flag
extern bool verboseMode;

// Logging functions
inline void log(const std::string& message) {
    if (verboseMode) {
        std::cout << "[LOG] " << message << std::endl;
    }
}

inline void logStep(const std::string& step, double value) {
    if (verboseMode) {
        std::cout << "[STEP] " << step << " = " << value << std::endl;
    }
}

inline void logStep(const std::string& step, const std::string& value) {
    if (verboseMode) {
        std::cout << "[STEP] " << step << ": " << value << std::endl;
    }
}

} // namespace Stats

#endif // LOGGER_H
