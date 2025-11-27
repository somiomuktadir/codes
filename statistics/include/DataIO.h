#ifndef DATAIO_H
#define DATAIO_H

#include <vector>
#include <string>
#include "Dataset.h"

namespace Stats {

class DataIO {
public:
    // CSV Import/Export
    static Dataset loadCSV(const std::string& filename, bool hasHeader = true, int column = 0);
    static std::vector<Dataset> loadMultiColumnCSV(const std::string& filename, bool hasHeader = true);
    
    static void saveCSV(const std::string& filename, const Dataset& data);
    static void saveCSV(const std::string& filename, const std::vector<Dataset>& datasets);
    static void saveCSV(const std::string& filename, const std::vector<double>& x, const std::vector<double>& y,
                       const std::string& xName = "X", const std::string& yName = "Y");
    
    // Data validation during import
    static bool validateFile(const std::string& filename);
    
private:
    static std::vector<std::string> splitLine(const std::string& line, char delimiter = ',');
    static bool isNumeric(const std::string& str);
};

} // namespace Stats

#endif // DATAIO_H
