#ifndef COMPLEX_MATRIX_H
#define COMPLEX_MATRIX_H

#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>

namespace LinAlg {

// Header-only complex matrix class for educational purposes
// For production use, consider integrating Eigen or another library
class ComplexMatrix {
private:
    std::vector<std::complex<double>> data;
    int rows;
    int cols;
    
public:
    // Constructors
    ComplexMatrix(int r, int c, std::complex<double> initialValue = std::complex<double>(0.0, 0.0))
        : rows(r), cols(c) {
        data.resize(r * c, initialValue);
    }
    
    // Accessors
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    
    std::complex<double>& at(int r, int c) {
        if (r < 0 || r >= rows || c < 0 || c >= cols) {
            throw std::out_of_range("ComplexMatrix index out of bounds");
        }
        return data[r * cols + c];
    }
    
    const std::complex<double>& at(int r, int c) const {
        if (r < 0 || r >= rows || c < 0 || c >= cols) {
            throw std::out_of_range("ComplexMatrix index out of bounds");
        }
        return data[r * cols + c];
    }
    
    std::complex<double>& operator()(int r, int c) {
        return data[r * cols + c];
    }
    
    const std::complex<double>& operator()(int r, int c) const {
        return data[r * cols + c];
    }
    
    // Arithmetic operations
    ComplexMatrix operator+(const ComplexMatrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for addition");
        }
        ComplexMatrix result(rows, cols);
        for (int i = 0; i < rows * cols; ++i) {
            result.data[i] = data[i] + other.data[i];
        }
        return result;
    }
    
    ComplexMatrix operator-(const ComplexMatrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for subtraction");
        }
        ComplexMatrix result(rows, cols);
        for (int i = 0; i < rows * cols; ++i) {
            result.data[i] = data[i] - other.data[i];
        }
        return result;
   }
    
    ComplexMatrix operator*(const ComplexMatrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Matrix dimensions incompatible for multiplication");
        }
        ComplexMatrix result(rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int k = 0; k < cols; ++k) {
                std::complex<double> aik = data[i * cols + k];
                for (int j = 0; j < other.cols; ++j) {
                    result.data[i * other.cols + j] += aik * other.data[k * other.cols + j];
                }
            }
        }
        return result;
    }
    
    ComplexMatrix operator*(std::complex<double> scalar) const {
        ComplexMatrix result(rows, cols);
        for (int i = 0; i < rows * cols; ++i) {
            result.data[i] = data[i] * scalar;
        }
        return result;
    }
    
    // Transpose
    ComplexMatrix transpose() const {
        ComplexMatrix result(cols, rows);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[j * rows + i] = data[i * cols + j];
            }
        }
        return result;
    }
    
    // Conjugate transpose (Hermitian transpose)
    ComplexMatrix conjugateTranspose() const {
        ComplexMatrix result(cols, rows);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.data[j * rows + i] = std::conj(data[i * cols + j]);
            }
        }
        return result;
    }
    
    // Check if Hermitian
    bool isHermitian(double tol = 1e-10) const {
        if (rows != cols) return false;
        for (int i = 0; i < rows; ++i) {
            for (int j = i; j < cols; ++j) {
                if (std::abs(data[i * cols + j] - std::conj(data[j * cols + i])) > tol) {
                    return false;
                }
            }
        }
        return true;
    }
    
    // Print
    void print() const {
        std::cout << std::fixed << std::setprecision(4);
        for (int i = 0; i < rows; ++i) {
            std::cout << "[ ";
            for (int j = 0; j < cols; ++j) {
                std::complex<double> val = data[i * cols + j];
                if (val.imag() >= 0) {
                    std::cout << val.real() << "+" << val.imag() << "i ";
                } else {
                    std::cout << val.real() << val.imag() << "i ";
                }
            }
            std::cout << "]" << std::endl;
        }
    }
    
    // Identity matrix
    static ComplexMatrix identity(int n) {
        ComplexMatrix result(n, n);
        for (int i = 0; i < n; ++i) {
            result.data[i * n + i] = std::complex<double>(1.0, 0.0);
        }
        return result;
    }
};

} // namespace LinAlg

#endif // COMPLEX_MATRIX_H
