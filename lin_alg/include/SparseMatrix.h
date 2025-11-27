#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "Matrix.h"
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace LinAlg {

// Compressed Sparse Row (CSR) format
class SparseMatrixCSR {
private:
    std::vector<double> values;      // Non-zero values
    std::vector<int> colIndices;     // Column indices of non-zero values
    std::vector<int> rowPtr;         // Row pointers (size = rows + 1)
    int rows;
    int cols;
    
public:
    // Constructors
    SparseMatrixCSR(int r, int c);
    
    // Convert from dense matrix
    static SparseMatrixCSR fromDense(const Matrix& dense, double threshold = 1e-10);
    
    // Convert to dense matrix
    Matrix toDense() const;
    
    // Accessors
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    int getNonZeros() const { return values.size(); }
    double getSparsity() const { 
        return 1.0 - (double)getNonZeros() / (double)(rows * cols); 
    }
    
    // Get element at (r, c) - O(log k) where k = non-zeros in row r
    double at(int r, int c) const;
    
    // Set element at (r, c) - use sparingly, prefer batch construction
    void set(int r, int c, double value);
    
    // Arithmetic operations
    SparseMatrixCSR operator+(const SparseMatrixCSR& other) const;
    SparseMatrixCSR operator*(double scalar) const;
    Matrix operator*(const Matrix& dense) const;  // Sparse * Dense -> Dense
    
    // Transpose to CSC
    class SparseMatrixCSC transpose() const;
    
    // Utilities
    void print() const;
    void printStats() const;
    
    // For iterative solvers - matrix-vector multiplication
    std::vector<double> multiply(const std::vector<double>& x) const;
    
    friend class SparseMatrixCSC;
};

// Compressed Sparse Column (CSC) format
class SparseMatrixCSC {
private:
    std::vector<double> values;      // Non-zero values
    std::vector<int> rowIndices;     // Row indices of non-zero values
    std::vector<int> colPtr;         // Column pointers (size = cols + 1)
    int rows;
    int cols;
    
public:
    // Constructors
    SparseMatrixCSC(int r, int c);
    
    // Convert from dense matrix
    static SparseMatrixCSC fromDense(const Matrix& dense, double threshold = 1e-10);
    
    // Convert to dense matrix
    Matrix toDense() const;
    
    // Accessors
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    int getNonZeros() const { return values.size(); }
    double getSparsity() const { 
        return 1.0 - (double)getNonZeros() / (double)(rows * cols); 
    }
    
    // Get element at (r, c) - O(log k) where k = non-zeros in column c
    double at(int r, int c) const;
    
    // Set element at (r, c)
    void set(int r, int c, double value);
    
    // Arithmetic operations
    SparseMatrixCSC operator+(const SparseMatrixCSC& other) const;
    SparseMatrixCSC operator*(double scalar) const;
    
    // Transpose to CSR
    SparseMatrixCSR transpose() const;
    
    // Utilities
    void print() const;
    void printStats() const;
    
    // For iterative solvers - matrix-vector multiplication
    std::vector<double> multiply(const std::vector<double>& x) const;
    
    friend class SparseMatrixCSR;
};

} // namespace LinAlg

#endif // SPARSE_MATRIX_H
