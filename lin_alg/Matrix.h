#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>

class Matrix {
private:
    std::vector<std::vector<double>> data;
    int rows;
    int cols;

public:
    // Constructors
    Matrix(int r, int c, double initialValue = 0.0);
    Matrix(const std::vector<std::vector<double>>& d);
    
    // Copy constructor and assignment
    Matrix(const Matrix& other) = default;
    Matrix& operator=(const Matrix& other) = default;
    
    // Move constructor and assignment for performance
    Matrix(Matrix&& other) noexcept = default;
    Matrix& operator=(Matrix&& other) noexcept = default;

    // Accessors
    int getRows() const;
    int getCols() const;
    double& at(int r, int c);
    const double& at(int r, int c) const;

    // Operations
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    Matrix transpose() const;
    
    // Analysis
    double trace() const;
    int rank() const;
    
    // Utilities
    static Matrix identity(int n);
    void print() const;
    
    // Friend functions for scalar multiplication (scalar * Matrix)
    friend Matrix operator*(double scalar, const Matrix& mat);
};

#endif // MATRIX_H
