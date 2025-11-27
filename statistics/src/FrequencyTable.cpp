#include "FrequencyTable.h"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace Stats {

void FrequencyTable::addInterval(double lower, double upper, int frequency) {
    FrequencyRow row;
    row.interval = {lower, upper};
    row.frequency = frequency;
    rows.push_back(row);
    updateCalculations();
}

int FrequencyTable::totalFrequency() const {
    int total = 0;
    for (const auto& row : rows) {
        total += row.frequency;
    }
    return total;
}

void FrequencyTable::updateCalculations() {
    int total = totalFrequency();
    int cumulative = 0;
    
    for (auto& row : rows) {
        cumulative += row.frequency;
        row.cumulativeFrequency = cumulative;
        row.relativeFrequency = (total > 0) ? static_cast<double>(row.frequency) / total : 0.0;
    }
}

FrequencyTable FrequencyTable::createFromData(const std::vector<double>& data, int numClasses) {
    FrequencyTable table;
    if (data.empty()) return table;
    
    auto minmax = std::minmax_element(data.begin(), data.end());
    double minVal = *minmax.first;
    double maxVal = *minmax.second;
    
    // Sturges' rule if numClasses not specified
    if (numClasses <= 0) {
        numClasses = static_cast<int>(1 + 3.322 * std::log10(data.size()));
        if (numClasses < 5) numClasses = 5; // Minimum 5 classes usually
    }
    
    double range = maxVal - minVal;
    // Add a small epsilon to range to ensure max value is included in last bin
    double classWidth = (range + 1e-9) / numClasses;
    
    // Create empty bins
    for (int i = 0; i < numClasses; ++i) {
        double lower = minVal + i * classWidth;
        double upper = minVal + (i + 1) * classWidth;
        table.addInterval(lower, upper, 0);
    }
    
    // Fill bins
    for (double val : data) {
        for (auto& row : table.rows) {
            if (row.interval.contains(val)) {
                row.frequency++;
                break;
            }
            // Handle the exact max value case (put in last bin)
            if (val == maxVal && row.interval.upper >= maxVal) {
                 // Check if it's the last bin
                 if (&row == &table.rows.back()) {
                     row.frequency++;
                     break;
                 }
            }
        }
    }
    
    table.updateCalculations();
    return table;
}

void FrequencyTable::print() const {
    std::cout << "\nFrequency Distribution Table" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;
    std::cout << std::setw(15) << "Class Interval" 
              << std::setw(10) << "Midpoint" 
              << std::setw(10) << "Freq" 
              << std::setw(12) << "Rel. Freq" 
              << std::setw(12) << "Cum. Freq" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;
    
    std::cout << std::fixed << std::setprecision(4);
    for (const auto& row : rows) {
        std::cout << "[" << std::setw(6) << row.interval.lower << " - " 
                  << std::setw(6) << row.interval.upper << ") "
                  << std::setw(10) << row.interval.midpoint()
                  << std::setw(10) << row.frequency
                  << std::setw(12) << row.relativeFrequency
                  << std::setw(12) << row.cumulativeFrequency << std::endl;
    }
    std::cout << "----------------------------------------------------------------" << std::endl;
    std::cout << "Total Frequency: " << totalFrequency() << std::endl;
}

} // namespace Stats
