#ifndef DATASET_H
#define DATASET_H

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>

namespace Stats {

class Dataset {
private:
    std::string name;
    std::vector<double> data;

public:
    // Constructors
    Dataset() : name("Unnamed") {}
    explicit Dataset(const std::string& n) : name(n) {}
    Dataset(const std::string& n, const std::vector<double>& d) : name(n), data(d) {}
    Dataset(const std::vector<double>& d) : name("Unnamed"), data(d) {}

    // Accessors
    const std::vector<double>& getData() const { return data; }
    std::vector<double>& getData() { return data; }
    const std::string& getName() const { return name; }
    void setName(const std::string& n) { name = n; }
    
    size_t size() const { return data.size(); }
    bool empty() const { return data.empty(); }
    
    double operator[](size_t i) const { return data[i]; }
    double& operator[](size_t i) { return data[i]; }
    
    // Data manipulation
    void add(double value) { data.push_back(value); }
    void clear() { data.clear(); }
    void sort() { std::sort(data.begin(), data.end()); }
    
    // Get sorted copy without modifying original
    std::vector<double> getSorted() const {
        std::vector<double> sorted = data;
        std::sort(sorted.begin(), sorted.end());
        return sorted;
    }
    
    // Data validation
    bool hasNaN() const {
        for (double val : data) {
            if (std::isnan(val)) return true;
        }
        return false;
    }
    
    bool hasInf() const {
        for (double val : data) {
            if (std::isinf(val)) return true;
        }
        return false;
    }
    
    // Remove outliers using IQR method
    void removeOutliers(double multiplier = 1.5);
    
    // Utilities
    void print() const;
    
    // Statistical summary (will use UnivariateStats)
    void printSummary() const;
};

} // namespace Stats

#endif // DATASET_H
