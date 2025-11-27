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
    Matrix transpose() const;
    
    // Element-wise operations
    Matrix hadamard(const Matrix& other) const;  // Element-wise multiplication
    Matrix applyFunction(std::function<double(double)> func) const;
    
    // Matrix manipulation
    void resize(int newRows, int newCols, double fillValue = 0.0);
    Matrix submatrix(int startRow, int startCol, int numRows, int numCols) const;
    static Matrix hstack(const Matrix& A, const Matrix& B);  // Horizontal concatenation
    static Matrix vstack(const Matrix& A, const Matrix& B);  // Vertical concatenation
    
    // File I/O
    void saveCSV(const std::string& filename) const;
    static Matrix loadCSV(const std::string& filename);
    
    // Analysis
    double trace() const;
    int rank() const;
    Matrix rref() const;
    double conditionNumber() const;
    
    // Utilities
    static Matrix identity(int n);
    void print() const;
    
    // Friend functions for scalar multiplication (scalar * Matrix)
    friend Matrix operator*(double scalar, const Matrix& mat);

    // Helpers
    static double epsilon() { return 1e-10; } // Placeholder for now, can be made dynamic later
};

} // namespace LinAlg

#endif // MATRIX_H
