#include "SparseMatrix.h"
#include "Logger.h"
#include <iomanip>

namespace LinAlg {

// ============================================================================
// SparseMatrixCSR Implementation
// ============================================================================

SparseMatrixCSR::SparseMatrixCSR(int r, int c) : rows(r), cols(c) {
    rowPtr.resize(r + 1, 0);
}

SparseMatrixCSR SparseMatrixCSR::fromDense(const Matrix& dense, double threshold) {
    Logger::getInstance().log("[STEP] Converting dense matrix to CSR format");
    
    SparseMatrixCSR sparse(dense.getRows(), dense.getCols());
    
    for (int i = 0; i < dense.getRows(); ++i) {
        sparse.rowPtr[i] = sparse.values.size();
        for (int j = 0; j < dense.getCols(); ++j) {
            double val = dense.at(i, j);
            if (std::abs(val) > threshold) {
                sparse.values.push_back(val);
                sparse.colIndices.push_back(j);
            }
        }
    }
    sparse.rowPtr[dense.getRows()] = sparse.values.size();
    
    Logger::getInstance().log("[STEP] CSR conversion complete. Sparsity: " + 
                              std::to_string(sparse.getSparsity() * 100) + "%");
    
    return sparse;
}

Matrix SparseMatrixCSR::toDense() const {
    Matrix dense(rows, cols, 0.0);
    
    for (int i = 0; i < rows; ++i) {
        for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
            dense.at(i, colIndices[k]) = values[k];
        }
    }
    
    return dense;
}

double SparseMatrixCSR::at(int r, int c) const {
    if (r < 0 || r >= rows || c < 0 || c >= cols) {
        throw std::out_of_range("CSR matrix index out of bounds");
    }
    
    // Binary search for column index in row r
    int start = rowPtr[r];
    int end = rowPtr[r + 1];
    
    auto it = std::lower_bound(colIndices.begin() + start, 
                               colIndices.begin() + end, c);
    
    if (it != colIndices.begin() + end && *it == c) {
        return values[it - colIndices.begin()];
    }
    
    return 0.0;
}

void SparseMatrixCSR::set(int r, int c, double value) {
    if (r < 0 || r >= rows || c < 0 || c >= cols) {
        throw std::out_of_range("CSR matrix index out of bounds");
    }
    
    // Find position in row r
    int start = rowPtr[r];
    int end = rowPtr[r + 1];
    
    auto it = std::lower_bound(colIndices.begin() + start, 
                               colIndices.begin() + end, c);
    int pos = it - colIndices.begin();
    
    if (it != colIndices.begin() + end && *it == c) {
        // Element exists, update it
        values[pos] = value;
    } else {
        // Insert new element
        values.insert(values.begin() + pos, value);
        colIndices.insert(colIndices.begin() + pos, c);
        
        // Update row pointers
        for (int i = r + 1; i <= rows; ++i) {
            rowPtr[i]++;
        }
    }
}

SparseMatrixCSR SparseMatrixCSR::operator+(const SparseMatrixCSR& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }
    
    SparseMatrixCSR result(rows, cols);
    
    for (int i = 0; i < rows; ++i) {
        result.rowPtr[i] = result.values.size();
        
        int k1 = rowPtr[i];
        int k2 = other.rowPtr[i];
        int end1 = rowPtr[i + 1];
        int end2 = other.rowPtr[i + 1];
        
        // Merge two sorted arrays
        while (k1 < end1 || k2 < end2) {
            if (k1 < end1 && (k2 >= end2 || colIndices[k1] < other.colIndices[k2])) {
                result.values.push_back(values[k1]);
                result.colIndices.push_back(colIndices[k1]);
                k1++;
            } else if (k2 < end2 && (k1 >= end1 || other.colIndices[k2] < colIndices[k1])) {
                result.values.push_back(other.values[k2]);
                result.colIndices.push_back(other.colIndices[k2]);
                k2++;
            } else {
                // Same column index
                double sum = values[k1] + other.values[k2];
                if (std::abs(sum) > 1e-10) {
                    result.values.push_back(sum);
                    result.colIndices.push_back(colIndices[k1]);
                }
                k1++;
                k2++;
            }
        }
    }
    result.rowPtr[rows] = result.values.size();
    
    return result;
}

SparseMatrixCSR SparseMatrixCSR::operator*(double scalar) const {
    SparseMatrixCSR result = *this;
    for (size_t i = 0; i < result.values.size(); ++i) {
        result.values[i] *= scalar;
    }
    return result;
}

Matrix SparseMatrixCSR::operator*(const Matrix& dense) const {
    if (cols != dense.getRows()) {
        throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
    }
    
    Matrix result(rows, dense.getCols(), 0.0);
    
    for (int i = 0; i < rows; ++i) {
        for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
            int j = colIndices[k];
            double val = values[k];
            
            for (int p = 0; p < dense.getCols(); ++p) {
                result.at(i, p) += val * dense.at(j, p);
            }
        }
    }
    
    return result;
}

std::vector<double> SparseMatrixCSR::multiply(const std::vector<double>& x) const {
    if ((int)x.size() != cols) {
        throw std::invalid_argument("Vector size must match matrix columns");
    }
    
    std::vector<double> result(rows, 0.0);
    
    for (int i = 0; i < rows; ++i) {
        for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
            result[i] += values[k] * x[colIndices[k]];
        }
    }
    
    return result;
}

void SparseMatrixCSR::print() const {
    std::cout << "CSR Matrix (" << rows << "x" << cols << ", " 
              << getNonZeros() << " non-zeros, "
              << std::fixed << std::setprecision(1) << getSparsity() * 100 << "% sparse)"
              << std::endl;
    
    // Print in dense format for small matrices
    if (rows <= 10 && cols <= 10) {
        Matrix dense = toDense();
        dense.print();
    } else {
        std::cout << "Matrix too large to display, showing non-zero entries:" << std::endl;
        for (int i = 0; i < rows; ++i) {
            for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                std::cout << "  (" << i << ", " << colIndices[k] << ") = " 
                          << std::fixed << std::setprecision(4) << values[k] << std::endl;
                if (k - rowPtr[i] >= 20) {
                    std::cout << "  ... (showing first 20 non-zeros)" << std::endl;
                    break;
                }
            }
        }
    }
}

void SparseMatrixCSR::printStats() const {
    std::cout << "CSR Matrix Statistics:" << std::endl;
    std::cout << "  Dimensions: " << rows << " x " << cols << std::endl;
    std::cout << "  Non-zeros: " << getNonZeros() << std::endl;
    std::cout << "  Sparsity: " << std::fixed << std::setprecision(2) 
              << getSparsity() * 100 << "%" << std::endl;
    std::cout << "  Memory (approx): " 
              << (values.size() * sizeof(double) + 
                  colIndices.size() * sizeof(int) + 
                  rowPtr.size() * sizeof(int)) / 1024.0 
              << " KB" << std::endl;
}

// ============================================================================
// SparseMatrixCSC Implementation
// ============================================================================

SparseMatrixCSC::SparseMatrixCSC(int r, int c) : rows(r), cols(c) {
    colPtr.resize(c + 1, 0);
}

SparseMatrixCSC SparseMatrixCSC::fromDense(const Matrix& dense, double threshold) {
    Logger::getInstance().log("[STEP] Converting dense matrix to CSC format");
    
    SparseMatrixCSC sparse(dense.getRows(), dense.getCols());
    
    for (int j = 0; j < dense.getCols(); ++j) {
        sparse.colPtr[j] = sparse.values.size();
        for (int i = 0; i < dense.getRows(); ++i) {
            double val = dense.at(i, j);
            if (std::abs(val) > threshold) {
                sparse.values.push_back(val);
                sparse.rowIndices.push_back(i);
            }
        }
    }
    sparse.colPtr[dense.getCols()] = sparse.values.size();
    
    Logger::getInstance().log("[STEP] CSC conversion complete. Sparsity: " + 
                              std::to_string(sparse.getSparsity() * 100) + "%");
    
    return sparse;
}

Matrix SparseMatrixCSC::toDense() const {
    Matrix dense(rows, cols, 0.0);
    
    for (int j = 0; j < cols; ++j) {
        for (int k = colPtr[j]; k < colPtr[j + 1]; ++k) {
            dense.at(rowIndices[k], j) = values[k];
        }
    }
    
    return dense;
}

double SparseMatrixCSC::at(int r, int c) const {
    if (r < 0 || r >= rows || c < 0 || c >= cols) {
        throw std::out_of_range("CSC matrix index out of bounds");
    }
    
    // Binary search for row index in column c
    int start = colPtr[c];
    int end = colPtr[c + 1];
    
    auto it = std::lower_bound(rowIndices.begin() + start, 
                               rowIndices.begin() + end, r);
    
    if (it != rowIndices.begin() + end && *it == r) {
        return values[it - rowIndices.begin()];
    }
    
    return 0.0;
}

void SparseMatrixCSC::set(int r, int c, double value) {
    if (r < 0 || r >= rows || c < 0 || c >= cols) {
        throw std::out_of_range("CSC matrix index out of bounds");
    }
    
    // Find position in column c
    int start = colPtr[c];
    int end = colPtr[c + 1];
    
    auto it = std::lower_bound(rowIndices.begin() + start, 
                               rowIndices.begin() + end, r);
    int pos = it - rowIndices.begin();
    
    if (it != rowIndices.begin() + end && *it == r) {
        // Element exists, update it
        values[pos] = value;
    } else {
        // Insert new element
        values.insert(values.begin() + pos, value);
        rowIndices.insert(rowIndices.begin() + pos, r);
        
        // Update column pointers
        for (int j = c + 1; j <= cols; ++j) {
            colPtr[j]++;
        }
    }
}

SparseMatrixCSC SparseMatrixCSC::operator+(const SparseMatrixCSC& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }
    
    SparseMatrixCSC result(rows, cols);
    
    for (int j = 0; j < cols; ++j) {
        result.colPtr[j] = result.values.size();
        
        int k1 = colPtr[j];
        int k2 = other.colPtr[j];
        int end1 = colPtr[j + 1];
        int end2 = other.colPtr[j + 1];
        
        // Merge two sorted arrays
        while (k1 < end1 || k2 < end2) {
            if (k1 < end1 && (k2 >= end2 || rowIndices[k1] < other.rowIndices[k2])) {
                result.values.push_back(values[k1]);
                result.rowIndices.push_back(rowIndices[k1]);
                k1++;
            } else if (k2 < end2 && (k1 >= end1 || other.rowIndices[k2] < rowIndices[k1])) {
                result.values.push_back(other.values[k2]);
                result.rowIndices.push_back(other.rowIndices[k2]);
                k2++;
            } else {
                // Same row index
                double sum = values[k1] + other.values[k2];
                if (std::abs(sum) > 1e-10) {
                    result.values.push_back(sum);
                    result.rowIndices.push_back(rowIndices[k1]);
                }
                k1++;
                k2++;
            }
        }
    }
    result.colPtr[cols] = result.values.size();
    
    return result;
}

SparseMatrixCSC SparseMatrixCSC::operator*(double scalar) const {
    SparseMatrixCSC result = *this;
    for (size_t i = 0; i < result.values.size(); ++i) {
        result.values[i] *= scalar;
    }
    return result;
}

std::vector<double> SparseMatrixCSC::multiply(const std::vector<double>& x) const {
    if ((int)x.size() != cols) {
        throw std::invalid_argument("Vector size must match matrix columns");
    }
    
    std::vector<double> result(rows, 0.0);
    
    for (int j = 0; j < cols; ++j) {
        for (int k = colPtr[j]; k < colPtr[j + 1]; ++k) {
            result[rowIndices[k]] += values[k] * x[j];
        }
    }
    
    return result;
}

void SparseMatrixCSC::print() const {
    std::cout << "CSC Matrix (" << rows << "x" << cols << ", " 
              << getNonZeros() << " non-zeros, "
              << std::fixed << std::setprecision(1) << getSparsity() * 100 << "% sparse)"
              << std::endl;
    
    // Print in dense format for small matrices
    if (rows <= 10 && cols <= 10) {
        Matrix dense = toDense();
        dense.print();
    } else {
        std::cout << "Matrix too large to display, showing non-zero entries:" << std::endl;
        for (int j = 0; j < cols; ++j) {
            for (int k = colPtr[j]; k < colPtr[j + 1]; ++k) {
                std::cout << "  (" << rowIndices[k] << ", " << j << ") = " 
                          << std::fixed << std::setprecision(4) << values[k] << std::endl;
                if (k - colPtr[j] >= 20) {
                    std::cout << "  ... (showing first 20 non-zeros)" << std::endl;
                    break;
                }
            }
        }
    }
}

void SparseMatrixCSC::printStats() const {
    std::cout << "CSC Matrix Statistics:" << std::endl;
    std::cout << "  Dimensions: " << rows << " x " << cols << std::endl;
    std::cout << "  Non-zeros: " << getNonZeros() << std::endl;
    std::cout << "  Sparsity: " << std::fixed << std::setprecision(2) 
              << getSparsity() * 100 << "%" << std::endl;
    std::cout << "  Memory (approx): " 
              << (values.size() * sizeof(double) + 
                  rowIndices.size() * sizeof(int) + 
                  colPtr.size() * sizeof(int)) / 1024.0 
              << " KB" << std::endl;
}

// Transpose conversions
SparseMatrixCSC SparseMatrixCSR::transpose() const {
    SparseMatrixCSC result(cols, rows);
    result.values = values;
    result.rowIndices = colIndices;
    result.colPtr = rowPtr;
    return result;
}

SparseMatrixCSR SparseMatrixCSC::transpose() const {
    SparseMatrixCSR result(cols, rows);
    result.values = values;
    result.colIndices = rowIndices;
    result.rowPtr = colPtr;
    return result;
}

} // namespace LinAlg
