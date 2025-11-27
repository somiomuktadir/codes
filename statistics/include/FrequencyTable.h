#ifndef FREQUENCY_TABLE_H
#define FREQUENCY_TABLE_H

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

namespace Stats {

struct Interval {
    double lower;
    double upper;
    
    double midpoint() const { return (lower + upper) / 2.0; }
    bool contains(double value) const { return value >= lower && value < upper; }
};

struct FrequencyRow {
    Interval interval;
    int frequency;
    double relativeFrequency;
    int cumulativeFrequency;
};

class FrequencyTable {
public:
    FrequencyTable() = default;
    
    // Add a class interval with frequency
    void addInterval(double lower, double upper, int frequency);
    
    // Generate from raw data
    static FrequencyTable createFromData(const std::vector<double>& data, int numClasses = -1);
    
    // Accessors
    const std::vector<FrequencyRow>& getRows() const { return rows; }
    size_t size() const { return rows.size(); }
    int totalFrequency() const;
    
    // Display
    void print() const;
    
private:
    std::vector<FrequencyRow> rows;
    void updateCalculations();
};

} // namespace Stats

#endif // FREQUENCY_TABLE_H
