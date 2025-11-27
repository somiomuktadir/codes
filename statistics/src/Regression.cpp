#include "Regression.h"
#include "UnivariateStats.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace Stats {

// Helper class for Matrix operations needed for regression
class Matrix {
public:
    std::vector<std::vector<double>> data;
    int rows, cols;

    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(r, std::vector<double>(c, 0.0));
    }

    static Matrix multiply(const Matrix& A, const Matrix& B) {
        if (A.cols != B.rows) throw std::invalid_argument("Matrix dimension mismatch");
        Matrix C(A.rows, B.cols);
        for (int i = 0; i < A.rows; ++i) {
            for (int j = 0; j < B.cols; ++j) {
                for (int k = 0; k < A.cols; ++k) {
                    C.data[i][j] += A.data[i][k] * B.data[k][j];
                }
            }
        }
        return C;
    }

    static Matrix transpose(const Matrix& A) {
        Matrix T(A.cols, A.rows);
        for (int i = 0; i < A.rows; ++i) {
            for (int j = 0; j < A.cols; ++j) {
                T.data[j][i] = A.data[i][j];
            }
        }
        return T;
    }

    // Gaussian elimination to solve Ax = b
    // Returns x vector
    static std::vector<double> solve(Matrix A, std::vector<double> b) {
        int n = A.rows;
        if (A.cols != n) throw std::invalid_argument("Matrix must be square for solving");
        if (b.size() != static_cast<size_t>(n)) throw std::invalid_argument("Vector size mismatch");

        // Augmented matrix [A|b]
        for (int i = 0; i < n; ++i) {
            A.data[i].push_back(b[i]);
        }

        // Forward elimination
        for (int i = 0; i < n; ++i) {
            // Pivot
            int maxRow = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(A.data[k][i]) > std::abs(A.data[maxRow][i])) {
                    maxRow = k;
                }
            }
            std::swap(A.data[i], A.data[maxRow]);

            if (std::abs(A.data[i][i]) < 1e-9) throw std::runtime_error("Matrix is singular");

            for (int k = i + 1; k < n; ++k) {
                double factor = A.data[k][i] / A.data[i][i];
                for (int j = i; j <= n; ++j) {
                    A.data[k][j] -= factor * A.data[i][j];
                }
            }
        }

        // Back substitution
        std::vector<double> x(n);
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0.0;
            for (int j = i + 1; j < n; ++j) {
                sum += A.data[i][j] * x[j];
            }
            x[i] = (A.data[i][n] - sum) / A.data[i][i];
        }

        return x;
    }
};

PolynomialRegressionResult Regression::polynomialRegression(
    const std::vector<double>& x, 
    const std::vector<double>& y, 
    int degree) {
    
    if (x.size() != y.size()) throw std::invalid_argument("Input sizes mismatch");
    if (x.size() <= static_cast<size_t>(degree)) throw std::invalid_argument("Not enough data points for degree");
    
    int n = x.size();
    int m = degree + 1; // Number of coefficients
    
    // Construct X matrix (Vandermonde)
    Matrix X(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            X.data[i][j] = std::pow(x[i], j);
        }
    }
    
    // Normal Equation: (X^T * X) * a = X^T * y
    Matrix XT = Matrix::transpose(X);
    Matrix XTX = Matrix::multiply(XT, X);
    
    // Compute X^T * y
    std::vector<double> XTy(m, 0.0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            XTy[i] += XT.data[i][j] * y[j];
        }
    }
    
    // Solve for coefficients
    std::vector<double> coeffs = Matrix::solve(XTX, XTy);
    
    // Calculate R-squared and Standard Error
    double ssTotal = 0.0;
    double ssRes = 0.0;
    double yMean = UnivariateStats::mean(y);
    
    for (int i = 0; i < n; ++i) {
        double pred = predictPolynomial(x[i], coeffs);
        ssTotal += (y[i] - yMean) * (y[i] - yMean);
        ssRes += (y[i] - pred) * (y[i] - pred);
    }
    
    PolynomialRegressionResult result;
    result.coefficients = coeffs;
    result.rSquared = 1.0 - (ssRes / ssTotal);
    result.standardError = std::sqrt(ssRes / (n - m));
    
    return result;
}

MultipleRegressionResult Regression::multipleLinearRegression(
    const std::vector<std::vector<double>>& x, 
    const std::vector<double>& y) {
    
    if (x.empty() || y.empty()) throw std::invalid_argument("Empty input");
    int n = y.size();
    int numFeatures = x.size(); // Number of independent variables
    
    if (x[0].size() != static_cast<size_t>(n)) throw std::invalid_argument("Feature dimension mismatch");
    
    int m = numFeatures + 1; // +1 for intercept
    
    // Construct X matrix (Design matrix)
    // Rows are samples, Cols are features (first col is 1s)
    Matrix X(n, m);
    for (int i = 0; i < n; ++i) {
        X.data[i][0] = 1.0; // Intercept term
        for (int j = 0; j < numFeatures; ++j) {
            X.data[i][j + 1] = x[j][i];
        }
    }
    
    // Normal Equation: (X^T * X) * b = X^T * y
    Matrix XT = Matrix::transpose(X);
    Matrix XTX = Matrix::multiply(XT, X);
    
    // Compute X^T * y
    std::vector<double> XTy(m, 0.0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            XTy[i] += XT.data[i][j] * y[j];
        }
    }
    
    // Solve for coefficients
    std::vector<double> coeffs = Matrix::solve(XTX, XTy);
    
    // Calculate R-squared and Standard Error
    double ssTotal = 0.0;
    double ssRes = 0.0;
    double yMean = UnivariateStats::mean(y);
    
    for (int i = 0; i < n; ++i) {
        // Prepare row for prediction
        std::vector<double> row(numFeatures);
        for(int j=0; j<numFeatures; ++j) row[j] = x[j][i];
        
        double pred = predictMultiple(row, coeffs);
        ssTotal += (y[i] - yMean) * (y[i] - yMean);
        ssRes += (y[i] - pred) * (y[i] - pred);
    }
    
    MultipleRegressionResult result;
    result.coefficients = coeffs;
    result.rSquared = 1.0 - (ssRes / ssTotal);
    result.standardError = std::sqrt(ssRes / (n - m));
    
    return result;
}

double Regression::predictPolynomial(double x, const std::vector<double>& coeffs) {
    double y = 0.0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        y += coeffs[i] * std::pow(x, i);
    }
    return y;
}

double Regression::predictMultiple(const std::vector<double>& x, const std::vector<double>& coeffs) {
    double y = coeffs[0]; // Intercept
    for (size_t i = 0; i < x.size(); ++i) {
        if (i + 1 < coeffs.size()) {
            y += coeffs[i + 1] * x[i];
        }
    }
    return y;
}

} // namespace Stats
