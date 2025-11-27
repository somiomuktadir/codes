#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>

namespace Utils {
    // Constants
    constexpr double EPSILON = 1e-9;
    constexpr double PI = 3.14159265358979323846;

    // Validation
    bool isInteger(double value);
    bool isPositive(double value);
    bool isNonNegative(double value);

    // Math helpers
    unsigned long long factorial(int n);
    unsigned long long nCr(int n, int r);
    unsigned long long nPr(int n, int r);

    // Advanced Math
    double gammaFunction(double x);
    double betaFunction(double x, double y);

    // Input helpers
    double getDoubleInput(const std::string& prompt);
    int getIntInput(const std::string& prompt);
    std::vector<double> getVectorInput(const std::string& prompt);
}

#endif // UTILS_H
