#include "LinearSolver.h"
#include "Logger.h"
#include "VectorOps.h"
#include <sstream>
#include <iomanip>

std::vector<double> LinearSolver::solve(const Matrix& A, const std::vector<double>& b) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    if (b.size() != n) throw std::invalid_argument("Vector dimension mismatch");

    // Create augmented matrix [A|b]
    // We'll work with a local copy of data for Gaussian elimination
    Matrix M = augment(A, b);
    
    // Gaussian Elimination
    Logger::getInstance().log("Starting Gaussian Elimination");
    Logger::getInstance().logStep("Initial Augmented Matrix:");
    // We can't easily print the matrix here without a helper, but we can log actions.

    for (int i = 0; i < n; ++i) {
        // Pivot
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(M.at(j, i)) > std::abs(M.at(pivot, i))) {
                pivot = j;
            }
        }
        
        // Swap rows
        if (pivot != i) {
            std::ostringstream oss;
            oss << "Swapping row " << i << " with row " << pivot;
            Logger::getInstance().log(oss.str());
            for (int k = 0; k <= n; ++k) {
                std::swap(M.at(i, k), M.at(pivot, k));
            }
        }
        
        if (std::abs(M.at(i, i)) < 1e-10) {
            throw std::runtime_error("Matrix is singular or nearly singular");
        }
        
        // Eliminate
        for (int j = i + 1; j < n; ++j) {
            double factor = M.at(j, i) / M.at(i, i);
            if (std::abs(factor) > 1e-10) {
                std::ostringstream oss;
                oss << "R" << j << " = R" << j << " - (" << factor << ") * R" << i;
                Logger::getInstance().log(oss.str());
                for (int k = i; k <= n; ++k) {
                    M.at(j, k) -= factor * M.at(i, k);
                }
            }
        }
    }
    
    // Back substitution
    Logger::getInstance().log("Starting Back Substitution");
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += M.at(i, j) * x[j];
        }
        x[i] = (M.at(i, n) - sum) / M.at(i, i);
    }
    
    return x;
}

double LinearSolver::determinant(const Matrix& A) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    
    Matrix M = A; // Copy
    double det = 1.0;
    Logger::getInstance().log("Calculating Determinant using Gaussian Elimination");
    
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(M.at(j, i)) > std::abs(M.at(pivot, i))) {
                pivot = j;
            }
        }
        
        if (pivot != i) {
            std::ostringstream oss;
            oss << "Swapping row " << i << " with row " << pivot << " (det sign flips)";
            Logger::getInstance().log(oss.str());
            for (int k = 0; k < n; ++k) {
                std::swap(M.at(i, k), M.at(pivot, k));
            }
            det *= -1;
        }
        
        if (std::abs(M.at(i, i)) < 1e-10) return 0.0;
        
        det *= M.at(i, i);
        
        for (int j = i + 1; j < n; ++j) {
            double factor = M.at(j, i) / M.at(i, i);
            if (std::abs(factor) > 1e-10) {
                 // Log only significant operations to avoid clutter
                 // Logger::getInstance().log("R" + std::to_string(j) + " -= " + std::to_string(factor) + " * R" + std::to_string(i));
            }
            for (int k = i; k < n; ++k) {
                M.at(j, k) -= factor * M.at(i, k);
            }
        }
    }
    
    return det;
}

Matrix LinearSolver::inverse(const Matrix& A) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    
    Matrix M = augmentIdentity(A);
    int cols = 2 * n;
    
    // Gauss-Jordan
    Logger::getInstance().log("Starting Gauss-Jordan for Inverse");
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(M.at(j, i)) > std::abs(M.at(pivot, i))) {
                pivot = j;
            }
        }
        
        if (pivot != i) {
            Logger::getInstance().log("Swapping row " + std::to_string(i) + " with row " + std::to_string(pivot));
            for (int k = 0; k < cols; ++k) {
                std::swap(M.at(i, k), M.at(pivot, k));
            }
        }
        
        double pivotVal = M.at(i, i);
        if (std::abs(pivotVal) < 1e-10) throw std::runtime_error("Matrix is singular");
        
        for (int k = 0; k < cols; ++k) {
            M.at(i, k) /= pivotVal;
        }
        Logger::getInstance().log("Dividing row " + std::to_string(i) + " by " + std::to_string(pivotVal));
        
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                double factor = M.at(j, i);
                for (int k = 0; k < cols; ++k) {
                    M.at(j, k) -= factor * M.at(i, k);
                }
            }
        }
    }
    
    // Extract inverse
    Matrix Inv(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Inv.at(i, j) = M.at(i, j + n);
        }
    }
    
    return Inv;
}

Matrix LinearSolver::power(const Matrix& A, int n) {
    if (A.getRows() != A.getCols()) throw std::invalid_argument("Matrix must be square");
    if (n < 0) return power(inverse(A), -n);
    if (n == 0) return Matrix::identity(A.getRows());
    
    Matrix res = Matrix::identity(A.getRows());
    Matrix base = A;
    while (n > 0) {
        if (n % 2 == 1) res = res * base;
        base = base * base;
        n /= 2;
    }
    return res;
}

Matrix LinearSolver::augment(const Matrix& A, const std::vector<double>& b) {
    int r = A.getRows();
    int c = A.getCols();
    Matrix M(r, c + 1);
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            M.at(i, j) = A.at(i, j);
        }
        M.at(i, c) = b[i];
    }
    return M;
}

Matrix LinearSolver::augmentIdentity(const Matrix& A) {
    int n = A.getRows();
    Matrix M(n, 2 * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            M.at(i, j) = A.at(i, j);
            M.at(i, j + n) = (i == j) ? 1.0 : 0.0;
        }
    }
    return M;
}

std::vector<double> LinearSolver::solveRefined(const Matrix& A, const std::vector<double>& b, int maxIter) {
    // Initial solution using standard Gaussian elimination
    std::vector<double> x = solve(A, b);
    
    Logger::getInstance().log("Starting iterative refinement");
    
    for (int iter = 0; iter < maxIter; ++iter) {
        // Compute residual: r = b - Ax
        std::vector<double> r(b.size());
        for (size_t i = 0; i < b.size(); ++i) {
            double sum = 0;
            for (int j = 0; j < A.getCols(); ++j) {
                sum += A.at(i, j) * x[j];
            }
            r[i] = b[i] - sum;
        }
        
        // Check if residual is small enough
        double residual_norm = VectorOps::norm(r);
        if (residual_norm < 1e-14) {
            Logger::getInstance().log("Converged at iteration " + std::to_string(iter) + ", residual: " + std::to_string(residual_norm));
            break;
        }
        
        // Solve A * delta_x = r
        std::vector<double> delta_x = solve(A, r);
        
        // Update solution: x = x + delta_x
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += delta_x[i];
        }
        
        Logger::getInstance().log("Iteration " + std::to_string(iter + 1) + ", residual norm: " + std::to_string(residual_norm));
    }
    
    return x;
}
