#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <string>

namespace LinAlg {

class Logger {
private:
    bool enabled;
    Logger() : enabled(false) {}
    
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
        log(message);
    }
};

} // namespace LinAlg

#endif // LOGGER_H
