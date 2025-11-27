#ifndef HYPOTHESIS_TESTING_H
#define HYPOTHESIS_TESTING_H

#include <vector>
#include <utility>
#include <string>

namespace Stats {

class HypothesisTesting {
public:
    // Test Result Structure
    struct TTestResult {
        double tStatistic;
        double pValue;
        double degreesOfFreedom;
        double confidenceIntervalLower;
        double confidenceIntervalUpper;
        bool significant;
        double meanDifference;
    };
    
    // One-Sample t-test
    // Tests if sample mean differs from a hypothesized value
    static TTestResult oneSampleTTest(const std::vector<double>& data, 
                                      double mu0, 
                                      double alpha = 0.05);
    
    // Two-Sample Independent t-test
    // Tests if two independent groups have different means
    static TTestResult twoSampleTTest(const std::vector<double>& group1,
                                      const std::vector<double>& group2,
                                      bool equalVariance = true,
                                      double alpha = 0.05);
    
    // Paired t-test
    // Tests if paired observations have different means
    static TTestResult pairedTTest(const std::vector<double>& before,
                                   const std::vector<double>& after,
                                   double alpha = 0.05);
    
    // Confidence Interval for Mean
    static std::pair<double, double> confidenceInterval(
        const std::vector<double>& data, 
        double confidence = 0.95);
    
    // Utilities
    static void printTTestResult(const TTestResult& result, const std::string& testName);
    
private:
    // t-distribution CDF approximation
    static double tDistributionCDF(double t, double df);
    
    // Get critical t-value for given alpha and df
    static double tCriticalValue(double alpha, double df, bool twoTailed = true);
    
    // Standard error calculations
    static double standardError(const std::vector<double>& data);
    static double pooledStandardError(const std::vector<double>& group1,
                                      const std::vector<double>& group2);
};

} // namespace Stats

#endif // HYPOTHESIS_TESTING_H
