#include "DataIO.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

namespace Stats {

Dataset DataIO::loadCSV(const std::string& filename, bool hasHeader, int column) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    Dataset dataset;
    std::string line;
    bool firstLine = true;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        // Skip header if present
        if (firstLine && hasHeader) {
            auto parts = splitLine(line);
            if (column < static_cast<int>(parts.size())) {
                dataset.setName(parts[column]);
            }
            firstLine = false;
            continue;
        }
        
        auto parts = splitLine(line);
        if (column >= 0 && column < static_cast<int>(parts.size())) {
            if (isNumeric(parts[column])) {
                dataset.add(std::stod(parts[column]));
            }
        }
        
        firstLine = false;
    }
    
    file.close();
    return dataset;
}

std::vector<Dataset> DataIO::loadMultiColumnCSV(const std::string& filename, bool hasHeader) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    std::vector<Dataset> datasets;
    std::string line;
    bool firstLine = true;
    std::vector<std::string> headers;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        auto parts = splitLine(line);
        
        if (firstLine) {
            if (hasHeader) {
                headers = parts;
                for (const auto& header : headers) {
                    datasets.push_back(Dataset(header));
                }
            } else {
                for (size_t i = 0; i < parts.size(); ++i) {
                    datasets.push_back(Dataset("Column" + std::to_string(i + 1)));
                    if (isNumeric(parts[i])) {
                        datasets[i].add(std::stod(parts[i]));
                    }
                }
            }
            firstLine = false;
            if (hasHeader) continue;
        }
        
        for (size_t i = 0; i < parts.size() && i < datasets.size(); ++i) {
            if (isNumeric(parts[i])) {
                datasets[i].add(std::stod(parts[i]));
            }
        }
    }
    
    file.close();
    return datasets;
}

void DataIO::saveCSV(const std::string& filename, const Dataset& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not create file: " + filename);
    }
    
    file << data.getName() << "\n";
    
    const auto& values = data.getData();
    for (double val : values) {
        file << val << "\n";
    }
    
    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

void DataIO::saveCSV(const std::string& filename, const std::vector<Dataset>& datasets) {
    if (datasets.empty()) {
        throw std::invalid_argument("No datasets to save");
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not create file: " + filename);
    }
    
    // Write headers
    for (size_t i = 0; i < datasets.size(); ++i) {
        file << datasets[i].getName();
        if (i < datasets.size() - 1) file << ",";
    }
    file << "\n";
    
    // Find max length
    size_t maxLen = 0;
    for (const auto& ds : datasets) {
        maxLen = std::max(maxLen, ds.size());
    }
    
    // Write data
    for (size_t row = 0; row < maxLen; ++row) {
        for (size_t col = 0; col < datasets.size(); ++col) {
            if (row < datasets[col].size()) {
                file << datasets[col][row];
            }
            if (col < datasets.size() - 1) file << ",";
        }
        file << "\n";
    }
    
    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

void DataIO::saveCSV(const std::string& filename, const std::vector<double>& x, const std::vector<double>& y,
                     const std::string& xName, const std::string& yName) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("X and Y must have the same size");
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not create file: " + filename);
    }
    
    file << xName << "," << yName << "\n";
    
    for (size_t i = 0; i < x.size(); ++i) {
        file << x[i] << "," << y[i] << "\n";
    }
    
    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

bool DataIO::validateFile(const std::string& filename) {
    std::ifstream file(filename);
    return file.is_open();
}

std::vector<std::string> DataIO::splitLine(const std::string& line, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string token;
    
    while (std::getline(ss, token, delimiter)) {
        // Trim whitespace
        token.erase(0, token.find_first_not_of(" \t\r\n"));
        token.erase(token.find_last_not_of(" \t\r\n") + 1);
        tokens.push_back(token);
    }
    
    return tokens;
}

bool DataIO::isNumeric(const std::string& str) {
    if (str.empty()) return false;
    
    char* end = nullptr;
    std::strtod(str.c_str(), &end);
    
    return end != str.c_str() && *end == '\0';
}

} // namespace Stats
