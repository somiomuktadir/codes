#include "HypothesisTesting.h"
#include "UnivariateStats.h"
#include "Logger.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace Stats {

// One-Sample t-test
HypothesisTesting::TTestResult HypothesisTesting::oneSampleTTest(
    const std::vector<double>& data, double mu0, double alpha) {
    
    if (data.empty()) throw std::invalid_argument("Cannot perform t-test on empty dataset");
    if (data.size() < 2) throw std::invalid_argument("Need at least 2 observations for t-test");
    
    log("Performing one-sample t-test...");
    
    TTestResult result;
    
    double mean = UnivariateStats::mean(data);
    double stdDev = UnivariateStats::standardDeviation(data, true);
    double se = stdDev / std::sqrt(data.size());
    
    result.tStatistic = (mean - mu0) / se;
    result.degreesOfFreedom = data.size() - 1;
    result.meanDifference = mean - mu0;
    
    // Two-tailed p-value
    double tAbs = std::abs(result.tStatistic);
    result.pValue = 2.0 * (1.0 - tDistributionCDF(tAbs, result.degreesOfFreedom));
    
    // Confidence interval
    double tCrit = tCriticalValue(alpha, result.degreesOfFreedom, true);
    result.confidenceIntervalLower = mean - tCrit * se;
    result.confidenceIntervalUpper = mean + tCrit * se;
    
    result.significant = result.pValue < alpha;
    
    logStep("Sample mean", mean);
    logStep("Hypothesized mean", mu0);
    logStep("t-statistic", result.tStatistic);
    logStep("p-value", result.pValue);
    
    return result;
}

// Two-Sample Independent t-test
HypothesisTesting::TTestResult HypothesisTesting::twoSampleTTest(
    const std::vector<double>& group1, const std::vector<double>& group2,
    bool equalVariance, double alpha) {
    
    if (group1.empty() || group2.empty()) {
        throw std::invalid_argument("Cannot perform t-test on empty groups");
    }
    if (group1.size() < 2 || group2.size() < 2) {
        throw std::invalid_argument("Each group needs at least 2 observations");
    }
    
    log(equalVariance ? "Performing two-sample t-test (equal variance)..." 
                      : "Performing Welch's t-test (unequal variance)...");
    
    TTestResult result;
    
    double mean1 = UnivariateStats::mean(group1);
    double mean2 = UnivariateStats::mean(group2);
    double var1 = UnivariateStats::variance(group1, true);
    double var2 = UnivariateStats::variance(group2, true);
    
    size_t n1 = group1.size();
    size_t n2 = group2.size();
    
    result.meanDifference = mean1 - mean2;
    
    if (equalVariance) {
        // Pooled variance
        double pooledVar = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2);
        double se = std::sqrt(pooledVar * (1.0/n1 + 1.0/n2));
        
        result.tStatistic = result.meanDifference / se;
        result.degreesOfFreedom = n1 + n2 - 2;
        
        // Confidence interval
        double tCrit = tCriticalValue(alpha, result.degreesOfFreedom, true);
        result.confidenceIntervalLower = result.meanDifference - tCrit * se;
        result.confidenceIntervalUpper = result.meanDifference + tCrit * se;
    } else {
        // Welch's t-test
        double se = std::sqrt(var1/n1 + var2/n2);
        
        result.tStatistic = result.meanDifference / se;
        
        // Welch-Satterthwaite degrees of freedom
        double num = std::pow(var1/n1 + var2/n2, 2);
        double denom = std::pow(var1/n1, 2)/(n1-1) + std::pow(var2/n2, 2)/(n2-1);
        result.degreesOfFreedom = num / denom;
        
        // Confidence interval
        double tCrit = tCriticalValue(alpha, result.degreesOfFreedom, true);
        result.confidenceIntervalLower = result.meanDifference - tCrit * se;
        result.confidenceIntervalUpper = result.meanDifference + tCrit * se;
    }
    
    // Two-tailed p-value
    double tAbs = std::abs(result.tStatistic);
    result.pValue = 2.0 * (1.0 - tDistributionCDF(tAbs, result.degreesOfFreedom));
    
    result.significant = result.pValue < alpha;
    
    logStep("Mean group 1", mean1);
    logStep("Mean group 2", mean2);
    logStep("Mean difference", result.meanDifference);
    logStep("t-statistic", result.tStatistic);
    logStep("p-value", result.pValue);
    
    return result;
}

// Paired t-test
HypothesisTesting::TTestResult HypothesisTesting::pairedTTest(
    const std::vector<double>& before, const std::vector<double>& after, double alpha) {
    
    if (before.size() != after.size()) {
        throw std::invalid_argument("Paired samples must have equal size");
    }
    if (before.empty()) {
        throw std::invalid_argument("Cannot perform t-test on empty datasets");
    }
    if (before.size() < 2) {
        throw std::invalid_argument("Need at least 2 pairs for paired t-test");
    }
    
    log("Performing paired t-test...");
    
    // Compute differences
    std::vector<double> differences;
    differences.reserve(before.size());
    
    for (size_t i = 0; i < before.size(); ++i) {
        differences.push_back(after[i] - before[i]);
    }
    
    // One-sample t-test on differences with mu0 = 0
    TTestResult result = oneSampleTTest(differences, 0.0, alpha);
    
    // Override mean difference to be meaningful
    double meanBefore = UnivariateStats::mean(before);
    double meanAfter = UnivariateStats::mean(after);
    result.meanDifference = meanAfter - meanBefore;
    
    logStep("Mean before", meanBefore);
    logStep("Mean after", meanAfter);
    
    return result;
}

// Confidence Interval for Mean
std::pair<double, double> HypothesisTesting::confidenceInterval(
    const std::vector<double>& data, double confidence) {
    
    if (data.empty()) throw std::invalid_argument("Cannot compute CI for empty dataset");
    if (data.size() < 2) throw std::invalid_argument("Need at least 2 observations for CI");
    if (confidence <= 0 || confidence >= 1) {
        throw std::invalid_argument("Confidence must be between 0 and 1");
    }
    
    log("Computing confidence interval...");
    
    double mean = UnivariateStats::mean(data);
    double stdDev = UnivariateStats::standardDeviation(data, true);
    double se = stdDev / std::sqrt(data.size());
    double df = data.size() - 1;
    
    double alpha = 1.0 - confidence;
    double tCrit = tCriticalValue(alpha, df, true);
    
    double lower = mean - tCrit * se;
    double upper = mean + tCrit * se;
    
    logStep("Confidence level", confidence * 100);
    logStep("CI lower", lower);
    logStep("CI upper", upper);
    
    return {lower, upper};
}

// Print t-test Result
void HypothesisTesting::printTTestResult(const TTestResult& result, const std::string& testName) {
    std::cout << "\n========== " << testName << " ==========" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "\nt-statistic:            " << result.tStatistic << std::endl;
    std::cout << "Degrees of freedom:     " << result.degreesOfFreedom << std::endl;
    std::cout << "p-value (two-tailed):   " << result.pValue << std::endl;
    std::cout << "Mean difference:        " << result.meanDifference << std::endl;
    
    std::cout << "\n95% Confidence Interval: [" 
              << result.confidenceIntervalLower << ", " 
              << result.confidenceIntervalUpper << "]" << std::endl;
    
    std::cout << "\nResult: " << (result.significant ? "SIGNIFICANT" : "NOT SIGNIFICANT") 
              << " (α = 0.05)" << std::endl;
    
    if (result.significant) {
        std::cout << "Decision: Reject the null hypothesis" << std::endl;
    } else {
        std::cout << "Decision: Fail to reject the null hypothesis" << std::endl;
    }
    
    std::cout << "========================================\n" << std::endl;
}

// t-distribution CDF (Approximation)
double HypothesisTesting::tDistributionCDF(double t, double df) {
    // Simple approximation using normal distribution for large df
    // For small df, use a polynomial approximation
    
    if (df > 30) {
        // Use normal approximation for large df
        return 0.5 * (1.0 + std::erf(t / std::sqrt(2.0)));
    }
    
    // Hill's approximation for t-distribution CDF
    double x = df / (df + t * t);
    double z = 0.5;
    
    // Incomplete beta function approximation
    // This is a simplified version - for production use, a proper implementation would be better
    if (t > 0) {
        z = 1.0 - 0.5 * std::pow(x, df/2.0);
    } else {
        z = 0.5 * std::pow(x, df/2.0);
    }
    
    return z;
}

// Critical t-value
double HypothesisTesting::tCriticalValue(double alpha, double df, bool twoTailed) {
    if (twoTailed) alpha = alpha / 2.0;
    
    // Approximation using inverse normal for large df
    if (df > 30) {
        // Normal approximation
        // For alpha = 0.025 (two-tailed 0.05), z ≈ 1.96
        return 1.96 + 0.5 / df;  // Simple correction for finite df
    }
    
    // Lookup table for common values (simplified)
    // In a production system, you'd use a more comprehensive table or algorithm
    if (df >= 20) return 2.086;
    if (df >= 15) return 2.131;
    if (df >= 10) return 2.228;
    if (df >= 5) return 2.571;
    
    return 2.776;  // Conservative estimate for small df
}

// Standard Error
double HypothesisTesting::standardError(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot compute SE for empty dataset");
    
    double stdDev = UnivariateStats::standardDeviation(data, true);
    return stdDev / std::sqrt(data.size());
}

// Pooled Standard Error
double HypothesisTesting::pooledStandardError(const std::vector<double>& group1,
                                               const std::vector<double>& group2) {
    if (group1.empty() || group2.empty()) {
        throw std::invalid_argument("Cannot compute pooled SE for empty groups");
    }
    
    double var1 = UnivariateStats::variance(group1, true);
    double var2 = UnivariateStats::variance(group2, true);
    size_t n1 = group1.size();
    size_t n2 = group2.size();
    
    double pooledVar = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2);
    return std::sqrt(pooledVar * (1.0/n1 + 1.0/n2));
}

} // namespace Stats
