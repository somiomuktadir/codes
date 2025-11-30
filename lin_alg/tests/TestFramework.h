#ifndef TEST_FRAMEWORK_H
#define TEST_FRAMEWORK_H

#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <cmath>
#include <sstream>

namespace TestFramework {

class TestRegistry {
public:
    static TestRegistry& getInstance() {
        static TestRegistry instance;
        return instance;
    }

    void registerTest(const std::string& name, std::function<void()> testFunc) {
        tests.push_back({name, testFunc});
    }

    int runAll() {
        int passed = 0;
        int failed = 0;
        std::cout << "Running " << tests.size() << " tests..." << std::endl;
        
        for (const auto& test : tests) {
            try {
                test.func();
                std::cout << "[PASS] " << test.name << std::endl;
                passed++;
            } catch (const std::exception& e) {
                std::cout << "[FAIL] " << test.name << ": " << e.what() << std::endl;
                failed++;
            } catch (...) {
                std::cout << "[FAIL] " << test.name << ": Unknown error" << std::endl;
                failed++;
            }
        }
        
        std::cout << "\nResults: " << passed << " passed, " << failed << " failed." << std::endl;
        return failed > 0 ? 1 : 0;
    }

private:
    struct Test {
        std::string name;
        std::function<void()> func;
    };
    std::vector<Test> tests;
};

class TestRegistrar {
public:
    TestRegistrar(const std::string& name, std::function<void()> func) {
        TestRegistry::getInstance().registerTest(name, func);
    }
};

#define TEST(name) \
    void name(); \
    TestFramework::TestRegistrar registrar_##name(#name, name); \
    void name()

// Assertion macros
#define ASSERT_TRUE(condition) \
    if (!(condition)) throw std::runtime_error("Assertion failed: " #condition);

#define ASSERT_FALSE(condition) \
    if (condition) throw std::runtime_error("Assertion failed: " #condition);

#define ASSERT_EQ(a, b) \
    if ((a) != (b)) { \
        std::stringstream ss; \
        ss << "Assertion failed: " << #a << " (" << (a) << ") != " << #b << " (" << (b) << ")"; \
        throw std::runtime_error(ss.str()); \
    }

#define ASSERT_NEAR(a, b, tol) \
    if (std::abs((a) - (b)) > (tol)) { \
        std::stringstream ss; \
        ss << "Assertion failed: |" << #a << " (" << (a) << ") - " << #b << " (" << (b) << ")| > " << tol; \
        throw std::runtime_error(ss.str()); \
    }

#define ASSERT_THROWS(statement, exceptionType) \
    try { \
        statement; \
        throw std::runtime_error("Expected exception " #exceptionType " not thrown"); \
    } catch (const exceptionType&) { \
        /* Expected */ \
    } catch (...) { \
        throw std::runtime_error("Caught unknown exception, expected " #exceptionType); \
    }

} // namespace TestFramework

#endif // TEST_FRAMEWORK_H
