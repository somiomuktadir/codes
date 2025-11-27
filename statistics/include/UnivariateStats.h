#ifndef UNIVARIATE_STATS_H
#define UNIVARIATE_STATS_H

#include <vector>
#include <map>
#include <string>
#include <cmath>

namespace Stats {

class UnivariateStats {
public:
    // Central Tendency Measures
    static double mean(const std::vector<double>& data);
    static double median(const std::vector<double>& data);
    static std::vector<double> mode(const std::vector<double>& data);
    static double geometricMean(const std::vector<double>& data);
    static double harmonicMean(const std::vector<double>& data);
    
    // Dispersion Measures
    static double variance(const std::vector<double>& data, bool sample = true);
    static double standardDeviation(const std::vector<double>& data, bool sample = true);
    static double range(const std::vector<double>& data);
    static double interquartileRange(const std::vector<double>& data);
    static double meanAbsoluteDeviation(const std::vector<double>& data);
    static double coefficientOfVariation(const std::vector<double>& data);
    
    // Distribution Shape
    static double skewness(const std::vector<double>& data);
    static double kurtosis(const std::vector<double>& data);
    
    // Position Measures
    static double quantile(const std::vector<double>& data, double q);
    static double percentile(const std::vector<double>& data, double p);
    static std::vector<double> quartiles(const std::vector<double>& data);
    
    // Frequency Analysis
    static std::map<double, int> frequencyDistribution(const std::vector<double>& data);
    static std::map<double, double> relativeFrequency(const std::vector<double>& data);
    static std::map<double, int> cumulativeFrequency(const std::vector<double>& data);
    
    // Summary Statistics
    struct FiveNumberSummary {
        double minimum;
        double q1;
        double median;
        double q3;
        double maximum;
    };
    static FiveNumberSummary fiveNumberSummary(const std::vector<double>& data);
    
    // Utility
    static void printSummary(const std::vector<double>& data);
    
private:
    // Helper for quantile calculation
    static double interpolateQuantile(const std::vector<double>& sorted, double position);
};

} // namespace Stats

#endif // UNIVARIATE_STATS_H
