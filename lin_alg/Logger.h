#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <string>
#include <vector>

class Logger {
public:
    static Logger& getInstance() {
        static Logger instance;
        return instance;
    }

    void enable() { enabled = true; }
    void disable() { enabled = false; }
    bool isEnabled() const { return enabled; }

    void log(const std::string& message) {
        if (enabled) {
            std::cout << "[STEP] " << message << std::endl;
        }
    }

    void logStep(const std::string& message) {
        if (enabled) {
            std::cout << "\n>>> " << message << std::endl;
        }
    }

private:
    Logger() : enabled(false) {}
    bool enabled;
};

#endif // LOGGER_H
