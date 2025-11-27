#ifndef GROUPED_DATASET_H
#define GROUPED_DATASET_H

#include "FrequencyTable.h"

namespace Stats {

class GroupedDataset {
public:
    GroupedDataset() = default;
    GroupedDataset(const FrequencyTable& table);
    
    void setTable(const FrequencyTable& table);
    const FrequencyTable& getTable() const { return table; }
    
    // Statistical Measures
    double mean() const;
    double median() const;
    std::vector<double> mode() const; // Can have multiple modes
    double variance(bool sample = true) const;
    double standardDeviation(bool sample = true) const;
    
private:
    FrequencyTable table;
};

} // namespace Stats

#endif // GROUPED_DATASET_H
