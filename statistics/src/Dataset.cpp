#include "Dataset.h"
#include "UnivariateStats.h"
#include <iostream>
#include <iomanip>

namespace Stats {

void Dataset::removeOutliers(double multiplier) {
    if (data.size() < 4) return; // Need at least 4 points for IQR method
    
    auto quartileVals = UnivariateStats::quartiles(data);
    double q1 = quartileVals[0];
    double q3 = quartileVals[2];
    double iqr = q3 - q1;
    
    double lowerBound = q1 - multiplier * iqr;
    double upperBound = q3 + multiplier * iqr;
    
    std::vector<double> filtered;
    for (double val : data) {
        if (val >= lowerBound && val <= upperBound) {
            filtered.push_back(val);
        }
    }
    
    data = filtered;
}

void Dataset::print() const {
    std::cout << "Dataset: " << name << " (" << data.size() << " values)" << std::endl;
    std::cout << "[ ";
    for (size_t i = 0; i < data.size() && i < 20; ++i) {
        std::cout << std::fixed << std::setprecision(4) << data[i];
        if (i < data.size() - 1 && i < 19) std::cout << ", ";
    }
    if (data.size() > 20) {
        std::cout << " ... (" << (data.size() - 20) << " more)";
    }
    std::cout << " ]" << std::endl;
}

void Dataset::printSummary() const {
    if (data.empty()) {
        std::cout << "Dataset is empty." << std::endl;
        return;
    }
    UnivariateStats::printSummary(data);
}

} // namespace Stats
