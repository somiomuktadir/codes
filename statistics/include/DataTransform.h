#ifndef DATA_TRANSFORM_H
#define DATA_TRANSFORM_H

#include <vector>
#include <string>

namespace Stats {

class DataTransform {
public:
    // Standardization (Z-score normalization)
    static std::vector<double> zScoreStandardize(const std::vector<double>& data);
    
    // Min-Max Normalization
    static std::vector<double> minMaxNormalize(const std::vector<double>& data,
                                                double newMin = 0.0, double newMax = 1.0);
    
    // Robust Statistics
    static double trimmedMean(const std::vector<double>& data, double trimProportion = 0.1);
    
    // Mathematical Transformations
    static std::vector<double> logTransform(const std::vector<double>& data);
    static std::vector<double> sqrtTransform(const std::vector<double>& data);
    static std::vector<double> powerTransform(const std::vector<double>& data, double power);
    
    // Reverse Transformations
    static std::vector<double> reverseZScore(const std::vector<double>& zScores,
                                              double originalMean, double originalStdDev);
    static std::vector<double> reverseMinMax(const std::vector<double>& normalized,
                                              double originalMin, double originalMax,
                                              double newMin = 0.0, double newMax = 1.0);
    
    // Utilities
    struct ScalingParams {
        double mean;
        double stdDev;
        double min;
        double max;
    };
    
    static ScalingParams getScalingParams(const std::vector<double>& data);
};

} // namespace Stats

#endif // DATA_TRANSFORM_H
