#include "Utils.h"
#include <iostream>
#include <sstream>
#include <algorithm>

namespace Utils {

    bool isInteger(double value) {
        return std::abs(value - std::round(value)) < EPSILON;
    }

    bool isPositive(double value) {
        return value > EPSILON;
    }

    bool isNonNegative(double value) {
        return value > -EPSILON;
    }

    unsigned long long factorial(int n) {
        if (n < 0) throw std::invalid_argument("Factorial is not defined for negative numbers.");
        if (n > 20) throw std::overflow_error("Factorial result too large for 64-bit integer.");
        
        unsigned long long result = 1;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }

    unsigned long long nCr(int n, int r) {
        if (r < 0 || r > n) return 0;
        if (r == 0 || r == n) return 1;
        if (r > n / 2) r = n - r;
        
        unsigned long long result = 1;
        for (int i = 1; i <= r; ++i) {
            result = result * (n - i + 1) / i;
        }
        return result;
    }

    unsigned long long nPr(int n, int r) {
        if (r < 0 || r > n) return 0;
        if (r == 0) return 1;
        
        unsigned long long result = 1;
        for (int i = 0; i < r; ++i) {
            result *= (n - i);
        }
        return result;
    }

    // Lanczos approximation for Gamma function
    double gammaFunction(double x) {
        if (x <= 0 && isInteger(x)) throw std::invalid_argument("Gamma function undefined for non-positive integers.");
        
        // For small x, use recurrence relation: gamma(x) = gamma(x+1)/x
        if (x < 0.5) {
            return Utils::PI / (std::sin(Utils::PI * x) * gammaFunction(1 - x));
        }

        static const double p[] = {
            0.99999999999980993, 676.5203681218851, -1259.1392167224028,
            771.32342877765313, -176.61502916214059, 12.507343278686905,
            -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
        };
        int g = 7;
        x -= 1;
        double a = p[0];
        double t = x + g + 0.5;
        for (int i = 1; i < 9; i++) {
            a += p[i] / (x + i);
        }
        return std::sqrt(2 * Utils::PI) * std::pow(t, x + 0.5) * std::exp(-t) * a;
    }

    double betaFunction(double x, double y) {
        if (x <= 0 || y <= 0) throw std::invalid_argument("Beta function defined for positive numbers.");
        return (gammaFunction(x) * gammaFunction(y)) / gammaFunction(x + y);
    }

    double getDoubleInput(const std::string& prompt) {
        double value;
        while (true) {
            std::cout << prompt;
            if (std::cin >> value) {
                return value;
            } else {
                std::cout << "Invalid input. Please enter a number." << std::endl;
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
    }

    int getIntInput(const std::string& prompt) {
        int value;
        while (true) {
            std::cout << prompt;
            if (std::cin >> value) {
                return value;
            } else {
                std::cout << "Invalid input. Please enter an integer." << std::endl;
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
    }

    std::vector<double> getVectorInput(const std::string& prompt) {
        std::cout << prompt << " (space-separated, end with non-number or newline): ";
        std::vector<double> values;
        double val;
        std::string line;
        std::getline(std::cin >> std::ws, line);
        std::stringstream ss(line);
        while (ss >> val) {
            values.push_back(val);
        }
        return values;
    }
}
