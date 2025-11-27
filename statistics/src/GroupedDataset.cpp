#include "GroupedDataset.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace Stats {

GroupedDataset::GroupedDataset(const FrequencyTable& table) : table(table) {}

void GroupedDataset::setTable(const FrequencyTable& table) {
    this->table = table;
}

double GroupedDataset::mean() const {
    double sumFx = 0.0;
    int N = table.totalFrequency();
    
    if (N == 0) return 0.0;
    
    for (const auto& row : table.getRows()) {
        sumFx += row.frequency * row.interval.midpoint();
    }
    
    return sumFx / N;
}

double GroupedDataset::median() const {
    int N = table.totalFrequency();
    if (N == 0) return 0.0;
    
    double halfN = N / 2.0;
    
    // Find median class
    const FrequencyRow* medianClass = nullptr;
    int prevCumFreq = 0;
    
    for (const auto& row : table.getRows()) {
        if (row.cumulativeFrequency >= halfN) {
            medianClass = &row;
            break;
        }
        prevCumFreq = row.cumulativeFrequency;
    }
    
    if (!medianClass) return 0.0;
    
    double L = medianClass->interval.lower;
    double f = medianClass->frequency;
    double h = medianClass->interval.upper - medianClass->interval.lower;
    double F = prevCumFreq;
    
    return L + ((halfN - F) / f) * h;
}

std::vector<double> GroupedDataset::mode() const {
    std::vector<double> modes;
    const auto& rows = table.getRows();
    if (rows.empty()) return modes;
    
    // Find max frequency
    int maxFreq = 0;
    for (const auto& row : rows) {
        if (row.frequency > maxFreq) maxFreq = row.frequency;
    }
    
    // Find all classes with max frequency
    for (size_t i = 0; i < rows.size(); ++i) {
        if (rows[i].frequency == maxFreq) {
            double L = rows[i].interval.lower;
            double h = rows[i].interval.upper - rows[i].interval.lower;
            double f1 = rows[i].frequency;
            double f0 = (i > 0) ? rows[i-1].frequency : 0.0;
            double f2 = (i < rows.size() - 1) ? rows[i+1].frequency : 0.0;
            
            double modeVal = L + ((f1 - f0) / (2*f1 - f0 - f2)) * h;
            modes.push_back(modeVal);
        }
    }
    
    return modes;
}

double GroupedDataset::variance(bool sample) const {
    double mu = mean();
    double sumFx2 = 0.0;
    int N = table.totalFrequency();
    
    if (N <= 1) return 0.0;
    
    for (const auto& row : table.getRows()) {
        double m = row.interval.midpoint();
        sumFx2 += row.frequency * (m - mu) * (m - mu);
    }
    
    return sumFx2 / (sample ? (N - 1) : N);
}

double GroupedDataset::standardDeviation(bool sample) const {
    return std::sqrt(variance(sample));
}

} // namespace Stats
