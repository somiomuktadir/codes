#include "UnivariateStats.h"
#include "Logger.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <set>

namespace Stats {

// Central Tendency Measures

double UnivariateStats::mean(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot compute mean of empty dataset");
    
    log("Computing mean...");
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double result = sum / data.size();
    
    logStep("Sum", sum);
    logStep("Count", data.size());
    logStep("Mean", result);
    
    return result;
}

double UnivariateStats::median(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot compute median of empty dataset");
    
    log("Computing median...");
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    
    size_t n = sorted.size();
    double result;
    
    if (n % 2 == 0) {
        result = (sorted[n/2 - 1] + sorted[n/2]) / 2.0;
        logStep("Even number of elements, averaging middle two", result);
    } else {
        result = sorted[n/2];
        logStep("Odd number of elements, taking middle value", result);
    }
    
    return result;
}

std::vector<double> UnivariateStats::mode(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot compute mode of empty dataset");
    
    log("Computing mode...");
    auto freq = frequencyDistribution(data);
    
    int maxFreq = 0;
    for (const auto& pair : freq) {
        if (pair.second > maxFreq) {
            maxFreq = pair.second;
        }
    }
    
    std::vector<double> modes;
    for (const auto& pair : freq) {
        if (pair.second == maxFreq) {
            modes.push_back(pair.first);
        }
    }
    
    logStep("Maximum frequency", maxFreq);
    logStep("Number of modes", modes.size());
    
    return modes;
}

double UnivariateStats::geometricMean(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot compute geometric mean of empty dataset");
    
    log("Computing geometric mean...");
    
    // Check for non-positive values
    for (double val : data) {
        if (val <= 0) {
            throw std::invalid_argument("Geometric mean requires all positive values");
        }
    }
    
    double logSum = 0.0;
    for (double val : data) {
        logSum += std::log(val);
    }
    
    double result = std::exp(logSum / data.size());
    logStep("Geometric mean", result);
    
    return result;
}

double UnivariateStats::harmonicMean(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot compute harmonic mean of empty dataset");
    
    log("Computing harmonic mean...");
    
    double reciprocalSum = 0.0;
    for (double val : data) {
        if (val == 0) {
            throw std::invalid_argument("Harmonic mean undefined for zero values");
        }
        reciprocalSum += 1.0 / val;
    }
    
    double result = data.size() / reciprocalSum;
    logStep("Harmonic mean", result);
    
    return result;
}

// Dispersion Measures

double UnivariateStats::variance(const std::vector<double>& data, bool sample) {
    if (data.empty()) throw std::invalid_argument("Cannot compute variance of empty dataset");
    if (sample && data.size() == 1) throw std::invalid_argument("Sample variance undefined for single value");
    
    log(sample ? "Computing sample variance (Welford's method)..." : "Computing population variance...");
    
    // Using Welford's online algorithm for numerical stability
    double m = 0.0;
    double s = 0.0;
    
    for (size_t i = 0; i < data.size(); ++i) {
        double delta = data[i] - m;
        m += delta / (i + 1);
        s += delta * (data[i] - m);
    }
    
    double divisor = sample ? (data.size() - 1) : data.size();
    double result = s / divisor;
    
    logStep("Variance", result);
    
    return result;
}

double UnivariateStats::standardDeviation(const std::vector<double>& data, bool sample) {
    log("Computing standard deviation...");
    double var = variance(data, sample);
    double result = std::sqrt(var);
    logStep("Standard deviation", result);
    return result;
}

double UnivariateStats::range(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot compute range of empty dataset");
    
    log("Computing range...");
    auto minmax = std::minmax_element(data.begin(), data.end());
    double result = *minmax.second - *minmax.first;
    
    logStep("Minimum", *minmax.first);
    logStep("Maximum", *minmax.second);
    logStep("Range", result);
    
    return result;
}

double UnivariateStats::interquartileRange(const std::vector<double>& data) {
    if (data.size() < 2) throw std::invalid_argument("Need at least 2 values for IQR");
    
    log("Computing interquartile range...");
    auto q = quartiles(data);
    double result = q[2] - q[0]; // Q3 - Q1
    
    logStep("Q1", q[0]);
    logStep("Q3", q[2]);
    logStep("IQR", result);
    
    return result;
}

double UnivariateStats::meanAbsoluteDeviation(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot compute MAD of empty dataset");
    
    log("Computing mean absolute deviation...");
    double m = mean(data);
    
    double sumAbsDev = 0.0;
    for (double val : data) {
        sumAbsDev += std::abs(val - m);
    }
    
    double result = sumAbsDev / data.size();
    logStep("Mean absolute deviation", result);
    
    return result;
}

double UnivariateStats::coefficientOfVariation(const std::vector<double>& data) {
    log("Computing coefficient of variation...");
    double m = mean(data);
    if (m == 0) throw std::invalid_argument("Coefficient of variation undefined when mean is zero");
    
    double sd = standardDeviation(data);
    double result = (sd / m) * 100.0;
    
    logStep("CV (%)", result);
    
    return result;
}

// Distribution Shape

double UnivariateStats::skewness(const std::vector<double>& data) {
    if (data.size() < 3) throw std::invalid_argument("Need at least 3 values for skewness");
    
    log("Computing skewness...");
    double m = mean(data);
    double sd = standardDeviation(data);
    
    if (sd == 0) return 0.0;
    
    double sum = 0.0;
    for (double val : data) {
        double z = (val - m) / sd;
        sum += z * z * z;
    }
    
    double result = sum / data.size();
    logStep("Skewness", result);
    
    return result;
}

double UnivariateStats::kurtosis(const std::vector<double>& data) {
    if (data.size() < 4) throw std::invalid_argument("Need at least 4 values for kurtosis");
    
    log("Computing excess kurtosis...");
    double m = mean(data);
    double sd = standardDeviation(data);
    
    if (sd == 0) return 0.0;
    
    double sum = 0.0;
    for (double val : data) {
        double z = (val - m) / sd;
        sum += z * z * z * z;
    }
    
    double result = (sum / data.size()) - 3.0; // Excess kurtosis
    logStep("Excess kurtosis", result);
    
    return result;
}

// Position Measures

double UnivariateStats::interpolateQuantile(const std::vector<double>& sorted, double position) {
    size_t n = sorted.size();
    
    if (position <= 0) return sorted[0];
    if (position >= n - 1) return sorted[n - 1];
    
    size_t lower = static_cast<size_t>(position);
    size_t upper = lower + 1;
    double fraction = position - lower;
    
    return sorted[lower] + fraction * (sorted[upper] - sorted[lower]);
}

double UnivariateStats::quantile(const std::vector<double>& data, double q) {
    if (data.empty()) throw std::invalid_argument("Cannot compute quantile of empty dataset");
    if (q < 0 || q > 1) throw std::invalid_argument("Quantile must be between 0 and 1");
    
    log("Computing quantile...");
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    
    double position = q * (sorted.size() - 1);
    double result = interpolateQuantile(sorted, position);
    
    logStep("Quantile position", position);
    logStep("Quantile value", result);
    
    return result;
}

double UnivariateStats::percentile(const std::vector<double>& data, double p) {
    return quantile(data, p / 100.0);
}

std::vector<double> UnivariateStats::quartiles(const std::vector<double>& data) {
    log("Computing quartiles...");
    return {quantile(data, 0.25), quantile(data, 0.5), quantile(data, 0.75)};
}

// Frequency Analysis

std::map<double, int> UnivariateStats::frequencyDistribution(const std::vector<double>& data) {
    std::map<double, int> freq;
    for (double val : data) {
        freq[val]++;
    }
    return freq;
}

std::map<double, double> UnivariateStats::relativeFrequency(const std::vector<double>& data) {
    if (data.empty()) return {};
    
    auto freq = frequencyDistribution(data);
    std::map<double, double> relFreq;
    
    for (const auto& pair : freq) {
        relFreq[pair.first] = static_cast<double>(pair.second) / data.size();
    }
    
    return relFreq;
}

std::map<double, int> UnivariateStats::cumulativeFrequency(const std::vector<double>& data) {
    auto freq = frequencyDistribution(data);
    std::map<double, int> cumFreq;
    
    int cumulative = 0;
    for (const auto& pair : freq) {
        cumulative += pair.second;
        cumFreq[pair.first] = cumulative;
    }
    
    return cumFreq;
}

// Summary Statistics

UnivariateStats::FiveNumberSummary UnivariateStats::fiveNumberSummary(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot compute summary of empty dataset");
    
    log("Computing five-number summary...");
    
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    
    FiveNumberSummary summary;
    summary.minimum = sorted[0];
    summary.maximum = sorted[sorted.size() - 1];
    summary.median = quantile(data, 0.5);
    summary.q1 = quantile(data, 0.25);
    summary.q3 = quantile(data, 0.75);
    
    return summary;
}

void UnivariateStats::printSummary(const std::vector<double>& data) {
    if (data.empty()) {
        std::cout << "No data to summarize." << std::endl;
        return;
    }
    
    std::cout << "\n========== Statistical Summary ==========" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    
    std::cout << "Count:                  " << data.size() << std::endl;
    
    auto fns = fiveNumberSummary(data);
    std::cout << "\n--- Five-Number Summary ---" << std::endl;
    std::cout << "Minimum:                " << fns.minimum << std::endl;
    std::cout << "Q1 (25th percentile):   " << fns.q1 << std::endl;
    std::cout << "Median (50th):          " << fns.median << std::endl;
    std::cout << "Q3 (75th percentile):   " << fns.q3 << std::endl;
    std::cout << "Maximum:                " << fns.maximum << std::endl;
    
    std::cout << "\n--- Central Tendency ---" << std::endl;
    std::cout << "Mean:                   " << mean(data) << std::endl;
    
    auto modes = mode(data);
    std::cout << "Mode:                   ";
    if (modes.size() == data.size()) {
        std::cout << "No mode (all unique)" << std::endl;
    } else {
        for (size_t i = 0; i < modes.size() && i < 5; ++i) {
            std::cout << modes[i];
            if (i < modes.size() - 1 && i < 4) std::cout << ", ";
        }
        if (modes.size() > 5) std::cout << "...";
        std::cout << std::endl;
    }
    
    std::cout << "\n--- Dispersion ---" << std::endl;
    std::cout << "Range:                  " << range(data) << std::endl;
    std::cout << "IQR:                    " << interquartileRange(data) << std::endl;
    std::cout << "Variance (sample):      " << variance(data, true) << std::endl;
    std::cout << "Std. Deviation (sample):" << standardDeviation(data, true) << std::endl;
    std::cout << "MAD:                    " << meanAbsoluteDeviation(data) << std::endl;
    
    try {
        std::cout << "CV:                     " << coefficientOfVariation(data) << "%" << std::endl;
    } catch (...) {
        std::cout << "CV:                     Undefined (mean is zero)" << std::endl;
    }
    
    if (data.size() >= 3) {
        std::cout << "\n--- Distribution Shape ---" << std::endl;
        std::cout << "Skewness:               " << skewness(data) << std::endl;
        if (data.size() >= 4) {
            std::cout << "Excess Kurtosis:        " << kurtosis(data) << std::endl;
        }
    }
    
    std::cout << "========================================\n" << std::endl;
}

} // namespace Stats
