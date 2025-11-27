#include "DataTransform.h"
#include "UnivariateStats.h"
#include "Logger.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Stats {

// Z-score Standardization
std::vector<double> DataTransform::zScoreStandardize(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot standardize empty dataset");
    
    log("Computing Z-score standardization...");
    
    double mean = UnivariateStats::mean(data);
    double stdDev = UnivariateStats::standardDeviation(data);
    
    if (stdDev == 0) {
        throw std::invalid_argument("Cannot standardize data with zero variance");
    }
    
    std::vector<double> zScores;
    zScores.reserve(data.size());
    
    for (double val : data) {
        zScores.push_back((val - mean) / stdDev);
    }
    
    logStep("Mean", mean);
    logStep("Std Dev", stdDev);
    log("Z-scores computed");
    
    return zScores;
}

// Min-Max Normalization
std::vector<double> DataTransform::minMaxNormalize(const std::vector<double>& data,
                                                    double newMin, double newMax) {
    if (data.empty()) throw std::invalid_argument("Cannot normalize empty dataset");
    
    log("Computing Min-Max normalization...");
    
    auto minmax = std::minmax_element(data.begin(), data.end());
    double oldMin = *minmax.first;
    double oldMax = *minmax.second;
    
    if (oldMin == oldMax) {
        throw std::invalid_argument("Cannot normalize data with no variation");
    }
    
    std::vector<double> normalized;
    normalized.reserve(data.size());
    
    double range = oldMax - oldMin;
    double newRange = newMax - newMin;
    
    for (double val : data) {
        double normalized_val = ((val - oldMin) / range) * newRange + newMin;
        normalized.push_back(normalized_val);
    }
    
    logStep("Original Min", oldMin);
    logStep("Original Max", oldMax);
    logStep("New Min", newMin);
    logStep("New Max", newMax);
    
    return normalized;
}

// Trimmed Mean
double DataTransform::trimmedMean(const std::vector<double>& data, double trimProportion) {
    if (data.empty()) throw std::invalid_argument("Cannot compute trimmed mean of empty dataset");
    if (trimProportion < 0 || trimProportion >= 0.5) {
        throw std::invalid_argument("Trim proportion must be in [0, 0.5)");
    }
    
    log("Computing trimmed mean...");
    
    std::vector<double> sorted = data;
    std::sort(sorted.begin(), sorted.end());
    
    size_t n = sorted.size();
    size_t trimCount = static_cast<size_t>(n * trimProportion);
    
    if (trimCount * 2 >= n) {
        throw std::invalid_argument("Too much trimming - no data left");
    }
    
    double sum = 0.0;
    size_t count = 0;
    
    for (size_t i = trimCount; i < n - trimCount; ++i) {
        sum += sorted[i];
        count++;
    }
    
    double result = sum / count;
    
    logStep("Trimmed count (each side)", static_cast<double>(trimCount));
    logStep("Remaining count", static_cast<double>(count));
    logStep("Trimmed mean", result);
    
    return result;
}

// Log Transform
std::vector<double> DataTransform::logTransform(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot transform empty dataset");
    
    log("Computing log transformation...");
    
    std::vector<double> transformed;
    transformed.reserve(data.size());
    
    for (double val : data) {
        if (val <= 0) {
            throw std::invalid_argument("Log transformation requires positive values");
        }
        transformed.push_back(std::log(val));
    }
    
    log("Log transformation complete");
    
    return transformed;
}

// Square Root Transform
std::vector<double> DataTransform::sqrtTransform(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot transform empty dataset");
    
    log("Computing square root transformation...");
    
    std::vector<double> transformed;
    transformed.reserve(data.size());
    
    for (double val : data) {
        if (val < 0) {
            throw std::invalid_argument("Square root transformation requires non-negative values");
        }
        transformed.push_back(std::sqrt(val));
    }
    
    log("Square root transformation complete");
    
    return transformed;
}

// Power Transform
std::vector<double> DataTransform::powerTransform(const std::vector<double>& data, double power) {
    if (data.empty()) throw std::invalid_argument("Cannot transform empty dataset");
    
    log("Computing power transformation...");
    
    std::vector<double> transformed;
    transformed.reserve(data.size());
    
    for (double val : data) {
        transformed.push_back(std::pow(val, power));
    }
    
    logStep("Power", power);
    log("Power transformation complete");
    
    return transformed;
}

// Reverse Z-score
std::vector<double> DataTransform::reverseZScore(const std::vector<double>& zScores,
                                                  double originalMean, double originalStdDev) {
    if (zScores.empty()) throw std::invalid_argument("Cannot reverse transform empty dataset");
    
    log("Reversing Z-score standardization...");
    
    std::vector<double> original;
    original.reserve(zScores.size());
    
    for (double z : zScores) {
        original.push_back(z * originalStdDev + originalMean);
    }
    
    return original;
}

// Reverse Min-Max
std::vector<double> DataTransform::reverseMinMax(const std::vector<double>& normalized,
                                                  double originalMin, double originalMax,
                                                  double newMin, double newMax) {
    if (normalized.empty()) throw std::invalid_argument("Cannot reverse transform empty dataset");
    
    log("Reversing Min-Max normalization...");
    
    std::vector<double> original;
    original.reserve(normalized.size());
    
    double originalRange = originalMax - originalMin;
    double newRange = newMax - newMin;
    
    for (double val : normalized) {
        double original_val = ((val - newMin) / newRange) * originalRange + originalMin;
        original.push_back(original_val);
    }
    
    return original;
}

// Get Scaling Parameters
DataTransform::ScalingParams DataTransform::getScalingParams(const std::vector<double>& data) {
    if (data.empty()) throw std::invalid_argument("Cannot get params from empty dataset");
    
    ScalingParams params;
    params.mean = UnivariateStats::mean(data);
    params.stdDev = UnivariateStats::standardDeviation(data);
    
    auto minmax = std::minmax_element(data.begin(), data.end());
    params.min = *minmax.first;
    params.max = *minmax.second;
    
    return params;
}

} // namespace Stats
