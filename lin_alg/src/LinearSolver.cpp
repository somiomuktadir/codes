#include "LinearSolver.h"
#include "Logger.h"
#include "VectorOps.h"
#include "Decomposer.h"
#include "Analysis.h"
#include <sstream>
#include <iomanip>

namespace LinAlg {

std::vector<double> LinearSolver::solve(const Matrix& A, const std::vector<double>& b) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    if (b.size() != static_cast<size_t>(n)) throw std::invalid_argument("Vector dimension mismatch");

    // Use PLU Decomposition: PA = LU => LUx = Pb
    auto [P, L, U] = Decomposer::PLU(A);
    
    // 1. Compute Pb
    std::vector<double> Pb = Analysis::transform(P, b);
    
    // 2. Solve Ly = Pb (Forward substitution)
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += L(i, j) * y[j];
        }
        y[i] = (Pb[i] - sum) / L(i, i);
    }
    
    // 3. Solve Ux = y (Back substitution)
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += U(i, j) * x[j];
        }
        if (std::abs(U(i, i)) < Matrix::epsilon()) {
            throw std::runtime_error("Matrix is singular");
        }
        x[i] = (y[i] - sum) / U(i, i);
    }
    
    return x;
}

double LinearSolver::determinant(const Matrix& A) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    
    // det(A) = det(P^-1 * L * U) = det(P)^-1 * det(L) * det(U)
    // det(P) is +1 or -1 depending on number of swaps.
    // det(L) is 1 (diagonal is all 1s).
    // det(U) is product of diagonal elements.
    
    auto [P, L, U] = Decomposer::PLU(A);
    
    double det = 1.0;
    for (int i = 0; i < n; ++i) {
        det *= U(i, i);
    }
    
    // Calculate sign of permutation
    // We can count swaps to find det(P). 
    // Since P is a permutation matrix, det(P) is (-1)^N_swaps.
    // We can find N_swaps by decomposing P.
    
    // Quick det(P) calculation:
    Matrix P_copy = P;
    int p_swaps = 0;
    for(int i=0; i<n; ++i) {
        if(P_copy(i,i) == 0) {
            for(int j=i+1; j<n; ++j) {
                if(P_copy(j,i) != 0) {
                    // Swap rows i and j
                    for(int k=0; k<n; ++k) std::swap(P_copy(i,k), P_copy(j,k));
                    p_swaps++;
                    break;
                }
            }
        }
    }
    
    if (p_swaps % 2 == 1) det *= -1;
    
    return det;
}

Matrix LinearSolver::inverse(const Matrix& A) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    
    // To find A^-1, we solve A * x_i = e_i for each column i of identity matrix
    Matrix Inv(n, n);
    Matrix I = Matrix::identity(n);
    
    // We can reuse the PLU decomposition for all columns
    auto [P, L, U] = Decomposer::PLU(A);
    
    for (int col = 0; col < n; ++col) {
        // e_i is the i-th column of I
        std::vector<double> e(n, 0.0);
        e[col] = 1.0;
        
        // Solve A * x = e
        // PA = LU => LUx = Pe
        
        // 1. Compute Pe
        std::vector<double> Pe = Analysis::transform(P, e);
        
        // 2. Solve Ly = Pe
        std::vector<double> y(n);
        for (int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < i; ++j) {
                sum += L(i, j) * y[j];
            }
            y[i] = (Pe[i] - sum) / L(i, i);
        }
        
        // 3. Solve Ux = y
        std::vector<double> x(n);
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0;
            for (int j = i + 1; j < n; ++j) {
                sum += U(i, j) * x[j];
            }
            if (std::abs(U(i, i)) < Matrix::epsilon()) {
                throw std::runtime_error("Matrix is singular");
            }
            x[i] = (y[i] - sum) / U(i, i);
        }
        
        // Store x in Inv
        for (int i = 0; i < n; ++i) {
            Inv(i, col) = x[i];
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
            M(i, j) = A(i, j);
        }
        M(i, c) = b[i];
    }
    return M;
}

Matrix LinearSolver::augmentIdentity(const Matrix& A) {
    int n = A.getRows();
    Matrix M(n, 2 * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            M(i, j) = A(i, j);
            M(i, j + n) = (i == j) ? 1.0 : 0.0;
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
                sum += A(i, j) * x[j];
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

std::vector<double> LinearSolver::leastSquares(const Matrix& A, const std::vector<double>& b) {
    // Solve Ax = b using QR decomposition
    // A = Q * R
    // Rx = Q^T * b
    
    int m = A.getRows();
    int n = A.getCols();
    
    if (m < n) throw std::invalid_argument("Underdetermined system (rows < cols) not supported by this least squares implementation");
    if (b.size() != static_cast<size_t>(m)) throw std::invalid_argument("Vector dimension mismatch");
    
    auto [Q, R] = Decomposer::QR(A);
    
    // Compute y = Q^T * b
    std::vector<double> y = Analysis::transform(Q.transpose(), b);
    
    // Solve Rx = y
    return solve(R, y);
}

std::vector<double> LinearSolver::solveCholesky(const Matrix& A, const std::vector<double>& b) {
    // For symmetric positive-definite systems: Ax = b
    // 1. Compute Cholesky decomposition: A = L * L^T
    // 2. Solve L * y = b (forward substitution)
    // 3. Solve L^T * x = y (back substitution)
    
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    if (b.size() != static_cast<size_t>(n)) throw std::invalid_argument("Vector dimension mismatch");
    
    Logger::getInstance().log("Solving using Cholesky decomposition");
    
    Matrix L = Decomposer::Cholesky(A);
    
    // Forward substitution: L * y = b
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += L(i, j) * y[j];
        }
        y[i] = (b[i] - sum) / L(i, i);
    }
    
    // Back substitution: L^T * x = y
    Matrix Lt = L.transpose();
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += Lt(i, j) * x[j];
        }
        x[i] = (y[i] - sum) / Lt(i, i);
    }
    
    return x;
}

} // namespace LinAlg
