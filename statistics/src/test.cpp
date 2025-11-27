#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "UnivariateStats.h"
#include "BivariateStats.h"

using namespace Stats;

bool approxEqual(double a, double b, double epsilon = 1e-6) {
    return std::abs(a - b) < epsilon;
}

int main() {
    std::cout << "Running Statistics Library Tests...\n" << std::endl;
    
    int passed = 0, total = 0;
    
    // Test 1: Univariate Mean
    {
        std::vector<double> data = {1, 2, 3, 4, 5};
        double result = UnivariateStats::mean(data);
        total++;
        if (approxEqual(result, 3.0)) {
            std::cout << "✓ Test 1 PASSED: Mean = " << result << std::endl;
            passed++;
        } else {
            std::cout << "✗ Test 1 FAILED: Expected 3.0, got " << result << std::endl;
        }
    }
    
    // Test 2: Univariate Median
    {
        std::vector<double> data = {1, 2, 3, 4, 5};
        double result = UnivariateStats::median(data);
        total++;
        if (approxEqual(result, 3.0)) {
            std::cout << "✓ Test 2 PASSED: Median = " << result << std::endl;
            passed++;
        } else {
            std::cout << "✗ Test 2 FAILED: Expected 3.0, got " << result << std::endl;
        }
    }
    
    // Test 3: Univariate Sample Variance
    {
        std::vector<double> data = {1, 2, 3, 4, 5};
        double result = UnivariateStats::variance(data, true);
        total++;
        // Variance = ((1-3)² + (2-3)² + (3-3)² + (4-3)² + (5-3)²) / 4 = 10/4 = 2.5
        if (approxEqual(result, 2.5)) {
            std::cout << "✓ Test 3 PASSED: Sample Variance = " << result << std::endl;
            passed++;
        } else {
            std::cout << "✗ Test 3 FAILED: Expected 2.5, got " << result << std::endl;
        }
    }
    
    // Test 4: Univariate Standard Deviation
    {
        std::vector<double> data = {1, 2, 3, 4, 5};
        double result = UnivariateStats::standardDeviation(data, true);
        total++;
        if (approxEqual(result, std::sqrt(2.5))) {
            std::cout << "✓ Test 4 PASSED: Std Dev = " << result << std::endl;
            passed++;
        } else {
            std::cout << "✗ Test 4 FAILED: Expected " << std::sqrt(2.5) << ", got " << result << std::endl;
        }
    }
    
    // Test 5: Univariate Range
    {
        std::vector<double> data = {1, 2, 3, 4, 5};
        double result = UnivariateStats::range(data);
        total++;
        if (approxEqual(result, 4.0)) {
            std::cout << "✓ Test 5 PASSED: Range = " << result << std::endl;
            passed++;
        } else {
            std::cout << "✗ Test 5 FAILED: Expected 4.0, got " << result << std::endl;
        }
    }
    
    // Test 6: Pearson Correlation (perfect positive)
    {
        std::vector<double> x = {1, 2, 3, 4, 5};
        std::vector<double> y = {2, 4, 6, 8, 10};
        double result = BivariateStats::pearsonCorrelation(x, y);
        total++;
        if (approxEqual(result, 1.0)) {
            std::cout << "✓ Test 6 PASSED: Pearson r (perfect) = " << result << std::endl;
            passed++;
        } else {
            std::cout << "✗ Test 6 FAILED: Expected 1.0, got " << result << std::endl;
        }
    }
    
    // Test 7: Covariance
    {
        std::vector<double> x = {1, 2, 3, 4, 5};
        std::vector<double> y = {2, 4, 5, 4, 6};
        double result = BivariateStats::covariance(x, y, true);
        total++;
        // Should be positive covariance
        if (result > 0) {
            std::cout << "✓ Test 7 PASSED: Covariance = " << result << std::endl;
            passed++;
        } else {
            std::cout << "✗ Test 7 FAILED: Expected positive covariance, got " << result << std::endl;
        }
    }
    
    // Test 8: Linear Regression
    {
        std::vector<double> x = {1, 2, 3, 4, 5};
        std::vector<double> y = {2, 4, 6, 8, 10};
        auto result = BivariateStats::linearRegression(x, y);
        total++;
        // For y = 2x, slope should be 2, intercept should be 0
        if (approxEqual(result.slope, 2.0) && approxEqual(result.intercept, 0.0)) {
            std::cout << "✓ Test 8 PASSED: Regression slope=" << result.slope 
                     << ", intercept=" << result.intercept << std::endl;
            passed++;
        } else {
            std::cout << "✗ Test 8 FAILED: Expected slope=2.0, intercept=0.0, got slope=" 
                     << result.slope << ", intercept=" << result.intercept << std::endl;
        }
    }
    
    // Test 9: R-squared (perfect fit)
    {
        std::vector<double> x = {1, 2, 3, 4, 5};
        std::vector<double> y = {2, 4, 6, 8, 10};
        auto result = BivariateStats::linearRegression(x, y);
        total++;
        if (approxEqual(result.rSquared, 1.0)) {
            std::cout << "✓ Test 9 PASSED: R-squared = " << result.rSquared << std::endl;
            passed++;
        } else {
            std::cout << "✗ Test 9 FAILED: Expected R² = 1.0, got " << result.rSquared << std::endl;
        }
    }
    
    // Test 10: Quartiles
    {
        std::vector<double> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        auto result = UnivariateStats::quartiles(data);
        total++;
        // Q1 should be around 3.25, Q2 around 5.5, Q3 around 7.75
        if (result.size() == 3 && result[0] > 2.0 && result[1] > 4.5 && result[2] > 7.0) {
            std::cout << "✓ Test 10 PASSED: Q1=" << result[0] 
                     << ", Q2=" << result[1] << ", Q3=" << result[2] << std::endl;
            passed++;
        } else {
            std::cout << "✗ Test 10 FAILED: Unexpected quartile values" << std::endl;
        }
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Results: " << passed << "/" << total << " tests passed" << std::endl;
    std::cout << "Success Rate: " << std::fixed << std::setprecision(1) 
             << (100.0 * passed / total) << "%" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    return (passed == total) ? 0 : 1;
}
