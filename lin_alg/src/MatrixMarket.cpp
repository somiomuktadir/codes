#include "MatrixMarket.h"
#include "Logger.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cctype>

namespace LinAlg {

void MatrixMarket::parseHeader(const std::string& line, Format& format, 
                                Field& field, Symmetry& symmetry) {
    std::istringstream iss(line);
    std::string token;
    
    // Expected format: %%MatrixMarket matrix <format> <field> <symmetry>
    iss >> token; // Should be "%%MatrixMarket"
    if (token != "%%MatrixMarket") {
        throw std::runtime_error("Invalid Matrix Market file: missing header");
    }
    
    iss >> token; // Should be "matrix"
    if (token != "matrix") {
        throw std::runtime_error("Only matrix type is supported");
    }
    
    // Parse format
    iss >> token;
    std::transform(token.begin(), token.end(), token.begin(), ::tolower);
    if (token == "coordinate") {
        format = Format::COORDINATE;
    } else if (token == "array") {
        format = Format::ARRAY;
    } else {
        throw std::runtime_error("Unknown format: " + token);
    }
    
    // Parse field
    iss >> token;
    std::transform(token.begin(), token.end(), token.begin(), ::tolower);
    if (token == "real") {
        field = Field::REAL;
    } else if (token == "integer") {
        field = Field::INTEGER;
    } else if (token == "complex") {
        field = Field::COMPLEX;
    } else if (token == "pattern") {
        field = Field::PATTERN;
    } else {
        throw std::runtime_error("Unknown field: " + token);
    }
    
    // Parse symmetry
    iss >> token;
    std::transform(token.begin(), token.end(), token.begin(), ::tolower);
    if (token == "general") {
        symmetry = Symmetry::GENERAL;
    } else if (token == "symmetric") {
        symmetry = Symmetry::SYMMETRIC;
    } else if (token == "skew-symmetric") {
        symmetry = Symmetry::SKEW_SYMMETRIC;
    } else if (token == "hermitian") {
        symmetry = Symmetry::HERMITIAN;
    } else {
        throw std::runtime_error("Unknown symmetry: " + token);
    }
}

void MatrixMarket::skipComments(std::ifstream& file) {
    std::streampos pos;
    std::string line;
    
    while (file) {
        pos = file.tellg();
        std::getline(file, line);
        if (!line.empty() && line[0] != '%') {
            file.seekg(pos);
            break;
        }
    }
}

void MatrixMarket::writeHeader(std::ofstream& file, Format format, 
                               Field field, Symmetry symmetry) {
    file << "%%MatrixMarket matrix ";
    
    if (format == Format::COORDINATE) {
        file << "coordinate ";
    } else {
        file << "array ";
    }
    
    if (field == Field::REAL) {
        file << "real ";
    } else if (field == Field::INTEGER) {
        file << "integer ";
    } else if (field == Field::COMPLEX) {
        file << "complex ";
    } else {
        file << "pattern ";
    }
    
    if (symmetry == Symmetry::GENERAL) {
        file << "general";
    } else if (symmetry == Symmetry::SYMMETRIC) {
        file << "symmetric";
    } else if (symmetry == Symmetry::SKEW_SYMMETRIC) {
        file << "skew-symmetric";
    } else {
        file << "hermitian";
    }
    
    file << std::endl;
}

Matrix MatrixMarket::loadDense(const std::string& filename) {
    Logger::getInstance().log("[STEP] Loading Matrix Market file: " + filename);
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string headerLine;
    std::getline(file, headerLine);
    
    Format format;
    Field field;
    Symmetry symmetry;
    parseHeader(headerLine, format, field, symmetry);
    
    if (field == Field::COMPLEX) {
        throw std::runtime_error("Complex matrices not yet supported in this function");
    }
    
    skipComments(file);
    
    int rows, cols;
    
    if (format == Format::ARRAY) {
        // Dense array format
        file >> rows >> cols;
        Logger::getInstance().log("[STEP] Array format: " + std::to_string(rows) + "x" + std::to_string(cols));
        
        Matrix mat(rows, cols);
        
        // Column-major order in Matrix Market
        for (int j = 0; j < cols; ++j) {
            for (int i = 0; i < rows; ++i) {
                double val;
                file >> val;
                mat.at(i, j) = val;
            }
        }
        
        return mat;
        
    } else {
        // Coordinate (sparse) format
        int nnz;
        file >> rows >> cols >> nnz;
        Logger::getInstance().log("[STEP] Coordinate format: " + std::to_string(rows) + "x" + 
                                  std::to_string(cols) + ", " + std::to_string(nnz) + " non-zeros");
        
        Matrix mat(rows, cols, 0.0);
        
        for (int k = 0; k < nnz; ++k) {
            int i, j;
            double val = 1.0;
            
            file >> i >> j;
            if (field != Field::PATTERN) {
                file >> val;
            }
            
            // Matrix Market uses 1-based indexing
            i--; j--;
            
            mat.at(i, j) = val;
            
            // Handle symmetry
            if (symmetry == Symmetry::SYMMETRIC && i != j) {
                mat.at(j, i) = val;
            } else if (symmetry == Symmetry::SKEW_SYMMETRIC && i != j) {
                mat.at(j, i) = -val;
            }
        }
        
        return mat;
    }
}

SparseMatrixCSR MatrixMarket::loadSparse(const std::string& filename) {
    Logger::getInstance().log("[STEP] Loading sparse Matrix Market file: " + filename);
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string headerLine;
    std::getline(file, headerLine);
    
    Format format;
    Field field;
    Symmetry symmetry;
    parseHeader(headerLine, format, field, symmetry);
    
    if (format != Format::COORDINATE) {
        throw std::runtime_error("Sparse loading requires coordinate format");
    }
    
    skipComments(file);
    
    int rows, cols, nnz;
    file >> rows >> cols >> nnz;
    
    // Build as dense first, then convert (simpler for handling symmetry)
    Matrix temp = loadDense(filename);
    return SparseMatrixCSR::fromDense(temp);
}

void MatrixMarket::saveDense(const std::string& filename, const Matrix& mat) {
    Logger::getInstance().log("[STEP] Saving matrix to Matrix Market file: " + filename);
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot create file: " + filename);
    }
    
    writeHeader(file, Format::ARRAY, Field::REAL, Symmetry::GENERAL);
    
    // Optional comment
    file << "% Generated by LinAlg Library" << std::endl;
    
    // Dimensions
    file << mat.getRows() << " " << mat.getCols() << std::endl;
    
    // Data in column-major order
    file << std::scientific << std::setprecision(15);
    for (int j = 0; j < mat.getCols(); ++j) {
        for (int i = 0; i < mat.getRows(); ++i) {
            file << mat.at(i, j) << std::endl;
        }
    }
    
    Logger::getInstance().log("[STEP] Matrix saved successfully");
}

void MatrixMarket::saveSparse(const std::string& filename, const SparseMatrixCSR& mat) {
    Logger::getInstance().log("[STEP] Saving sparse matrix to Matrix Market file: " + filename);
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot create file: " + filename);
    }
    
    writeHeader(file, Format::COORDINATE, Field::REAL, Symmetry::GENERAL);
    
    // Optional comment
    file << "% Generated by LinAlg Library (Sparse CSR)" << std::endl;
    
    // Dimensions and non-zeros
    file << mat.getRows() << " " << mat.getCols() << " " << mat.getNonZeros() << std::endl;
    
    // Convert to dense temporarily to write (inefficient but simple)
    Matrix dense = mat.toDense();
    
    file << std::scientific << std::setprecision(15);
    for (int i = 0; i < dense.getRows(); ++i) {
        for (int j = 0; j < dense.getCols(); ++j) {
            double val = dense.at(i, j);
            if (std::abs(val) > 1e-10) {
                // Matrix Market uses 1-based indexing
                file << (i + 1) << " " << (j + 1) << " " << val << std::endl;
            }
        }
    }
    
    Logger::getInstance().log("[STEP] Sparse matrix saved successfully");
}

MatrixMarket::Format MatrixMarket::getFormat(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string headerLine;
    std::getline(file, headerLine);
    
    Format format;
    Field field;
    Symmetry symmetry;
    parseHeader(headerLine, format, field, symmetry);
    
    return format;
}

} // namespace LinAlg
