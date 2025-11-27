#include "TimeSeries.h"
#include "UnivariateStats.h"
#include <cmath>
#include <stdexcept>

namespace Stats {

double TimeSeries::autocorrelation(const std::vector<double>& data, int lag) {
    if (data.empty()) throw std::invalid_argument("Empty dataset");
    if (lag < 0) throw std::invalid_argument("Lag must be non-negative");
    if (static_cast<size_t>(lag) >= data.size()) return 0.0; // Not enough data
    
    double mean = UnivariateStats::mean(data);
    double numerator = 0.0;
    double denominator = 0.0;
    
    size_t n = data.size();
    
    for (size_t t = 0; t < n; ++t) {
        denominator += (data[t] - mean) * (data[t] - mean);
    }
    
    if (denominator == 0) return 0.0; // Constant signal
    
    for (size_t t = 0; t < n - lag; ++t) {
        numerator += (data[t] - mean) * (data[t + lag] - mean);
    }
    
    return numerator / denominator;
}

std::vector<double> TimeSeries::calculateACF(const std::vector<double>& data, int maxLag) {
    std::vector<double> acf;
    if (data.empty()) return acf;
    
    if (maxLag < 0) maxLag = 0;
    if (static_cast<size_t>(maxLag) >= data.size()) maxLag = data.size() - 1;
    
    for (int k = 0; k <= maxLag; ++k) {
        acf.push_back(autocorrelation(data, k));
    }
    
    return acf;
}

} // namespace Stats
