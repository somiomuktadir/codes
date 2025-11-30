#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <initializer_list>
#include <functional>
#include <string>

namespace LinAlg {

class Matrix {
private:
    std::vector<double> data;  // Flat storage for better cache locality
    int rows;
    int cols;

public:
    // Constructors
    Matrix(int r, int c, double initialValue = 0.0);
    Matrix(const std::vector<std::vector<double>>& d);
    Matrix(std::initializer_list<std::initializer_list<double>> list);
    
    // Copy constructor and assignment
    Matrix(const Matrix& other) = default;
    Matrix& operator=(const Matrix& other) = default;
    
    // Move constructor and assignment for performance
    Matrix(Matrix&& other) noexcept = default;
    Matrix& operator=(Matrix&& other) noexcept = default;

    // Accessors
    int getRows() const;
    int getCols() const;
    
    // Checked access (bounds checking)
    double& at(int r, int c);
    const double& at(int r, int c) const;
    
    // Unchecked access (faster, use in tight loops)
    double& operator()(int r, int c);
    const double& operator()(int r, int c) const;

    // Comparison operators
    bool operator==(const Matrix& other) const;
    bool operator!=(const Matrix& other) const;

    // Arithmetic operations
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    
    // Compound assignment operators (New)
    Matrix& operator+=(const Matrix& other);
    Matrix& operator-=(const Matrix& other);
    Matrix& operator*=(double scalar);
    // Note: operator*=(const Matrix&) is not typically element-wise, so we avoid it to prevent confusion, 
    // or implement it as matrix multiplication (which changes dimensions, so maybe not safe for inplace).

    Matrix transpose() const;
    
    // Element-wise operations
    Matrix hadamard(const Matrix& other) const;  // Element-wise multiplication
    Matrix applyFunction(std::function<double(double)> func) const;
    
    // Matrix manipulation
    void resize(int newRows, int newCols, double fillValue = 0.0);
    Matrix submatrix(int startRow, int startCol, int numRows, int numCols) const;
    static Matrix hstack(const Matrix& A, const Matrix& B);  // Horizontal concatenation
    static Matrix vstack(const Matrix& A, const Matrix& B);  // Vertical concatenation
    
    // Matrix Norms
    double frobeniusNorm() const;
    double l1Norm() const;   // Maximum absolute column sum
    double lInfNorm() const; // Maximum absolute row sum
    
    // Advanced Operations
    static Matrix kroneckerProduct(const Matrix& A, const Matrix& B);
    
    // Analysis
    double trace() const;
    int rank() const;
    Matrix rref() const;
    double conditionNumber() const;
    
    // Helpers (New)
    bool isSquare() const;
    bool isSymmetric(double tol = 1e-10) const;
    
    // Iterators (New)
    std::vector<double>::iterator begin() { return data.begin(); }
    std::vector<double>::iterator end() { return data.end(); }
    std::vector<double>::const_iterator begin() const { return data.begin(); }
    std::vector<double>::const_iterator end() const { return data.end(); }
    
    // Utilities
    static Matrix identity(int n);
    void print() const;
    
    // Friend functions for scalar multiplication (scalar * Matrix)
    friend Matrix operator*(double scalar, const Matrix& mat);

    // Helpers
    static double epsilon() { return 1e-10; }
};

} // namespace LinAlg

#endif // MATRIX_H
