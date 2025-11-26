#include "Decomposer.h"
#include "Decomposer.h"
#include "LinearSolver.h" // For inverse if needed, or just math
#include "Logger.h"
#include <algorithm>
#include <sstream>
#include <iomanip>

std::pair<Matrix, Matrix> Decomposer::LU(const Matrix& A) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");

    Matrix L = Matrix::identity(n);
    Matrix U = A;
    Logger::getInstance().log("Starting LU Decomposition");

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(U.at(i, i)) < 1e-10) throw std::runtime_error("Pivot is zero");
            double factor = U.at(j, i) / U.at(i, i);
            L.at(j, i) = factor;
            // Logger::getInstance().log("L(" + std::to_string(j) + "," + std::to_string(i) + ") = " + std::to_string(factor));
            for (int k = i; k < n; ++k) {
                U.at(j, k) -= factor * U.at(i, k);
            }
        }
    }
    Logger::getInstance().log("LU Decomposition completed");
    return std::make_pair(L, U);
}

Matrix Decomposer::Cholesky(const Matrix& A) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    
    Matrix L(n, n);
    Logger::getInstance().log("Starting Cholesky Decomposition");
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += L.at(i, k) * L.at(j, k);
            }
            
            if (i == j) {
                double val = A.at(i, i) - sum;
                if (val <= 0) throw std::runtime_error("Matrix is not positive definite");
                L.at(i, j) = std::sqrt(val);
            } else {
                L.at(i, j) = (A.at(i, j) - sum) / L.at(j, j);
            }
        }
    }
    return L;
}

std::pair<Matrix, Matrix> Decomposer::QR(const Matrix& A) {
    int m = A.getRows();
    int n = A.getCols();
    
    Matrix Q(m, n);
    Matrix R(n, n);
    Logger::getInstance().log("Starting QR Decomposition (Modified Gram-Schmidt)");
    
    // Modified Gram-Schmidt is more numerically stable than Classical GS
    std::vector<std::vector<double>> qCols(n, std::vector<double>(m));
    
    // Initialize Q columns with A columns
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            qCols[j][i] = A.at(i, j);
        }
    }
    
    // Modified Gram-Schmidt process
    for (int j = 0; j < n; ++j) {
        // Orthogonalize against all previous vectors
        for (int i = 0; i < j; ++i) {
            // r_ij = <q_i, a_j>
            double r_ij = VectorOps::dot(qCols[i], qCols[j]);
            R.at(i, j) = r_ij;
            
            // a_j = a_j - r_ij * q_i
            for (int k = 0; k < m; ++k) {
                qCols[j][k] -= r_ij * qCols[i][k];
            }
        }
        
        // Normalize q_j
        double r_jj = VectorOps::norm(qCols[j]);
        R.at(j, j) = r_jj;
        
        if (std::abs(r_jj) > 1e-10) {
            for (int k = 0; k < m; ++k) {
                qCols[j][k] /= r_jj;
            }
        } else {
            // Handle linearly dependent vectors
            Logger::getInstance().log("Warning: Linearly dependent column detected at index " + std::to_string(j));
        }
    }
    
    // Copy to Q matrix
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            Q.at(i, j) = qCols[j][i];
        }
    }
    
    return std::make_pair(Q, R);
}

std::pair<Matrix, Matrix> Decomposer::Eigen(const Matrix& A, int maxIter, double tol) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    
    Matrix T = A;
    Matrix Q_total = Matrix::identity(n);
    Logger::getInstance().log("Starting Eigen Decomposition (QR Algorithm)");
    
    for (int iter = 0; iter < maxIter; ++iter) {
        auto [Q, R] = QR(T);
        T = R * Q;
        Q_total = Q_total * Q;
        
        // Check convergence (off-diagonal elements close to 0)
        double offDiagonalSum = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                offDiagonalSum += std::abs(T.at(i, j));
            }
        }
        if (offDiagonalSum < tol) {
            Logger::getInstance().log("Converged at iteration " + std::to_string(iter));
            break;
        }
    }
    
    // Eigenvalues are on diagonal of T
    Matrix Eigenvalues(n, n); // Diagonal matrix
    for(int i=0; i<n; ++i) Eigenvalues.at(i,i) = T.at(i,i);

    // For symmetric matrices, Q_total contains eigenvectors.
    // For non-symmetric, we need to solve (T - lambda*I)x = 0 to get eigenvectors of T,
    // then transform back with Q_total * x.
    
    // Check if symmetric
    bool symmetric = true;
    for(int i=0; i<n; ++i) {
        for(int j=i+1; j<n; ++j) {
            if(std::abs(A.at(i,j) - A.at(j,i)) > 1e-5) {
                symmetric = false;
                break;
            }
        }
    }
    
    if (symmetric) {
        Logger::getInstance().log("Matrix is symmetric, Q contains eigenvectors");
        return {Eigenvalues, Q_total};
    }
    
    Logger::getInstance().log("Matrix is non-symmetric, computing eigenvectors via back-substitution");
    Matrix Eigenvectors(n, n);
    
    for (int k = 0; k < n; ++k) {
        double lambda = T.at(k, k);
        // Solve (T - lambda*I)x = 0
        // Since T is upper triangular, we can use back-substitution.
        // We want a non-trivial solution.
        // Let x[k] = 1 (arbitrary scale), then solve for x[0]...x[k-1].
        // x[k+1]...x[n-1] are 0.
        
        std::vector<double> x(n, 0.0);
        x[k] = 1.0;
        
        for (int i = k - 1; i >= 0; --i) {
            double sum = 0;
            for (int j = i + 1; j <= k; ++j) {
                sum += T.at(i, j) * x[j];
            }
            // (T[i][i] - lambda) * x[i] + sum = 0
            // x[i] = -sum / (T[i][i] - lambda)
            double denom = T.at(i, i) - lambda;
            if (std::abs(denom) < 1e-10) {
                // If we hit a zero divisor, it means repeated eigenvalue or similar.
                // This simple implementation might fail for defective matrices.
                x[i] = 0; // Fallback
            } else {
                x[i] = -sum / denom;
            }
        }
        
        // Transform back: v = Q * x
        for (int i = 0; i < n; ++i) {
            double val = 0;
            for (int j = 0; j < n; ++j) {
                val += Q_total.at(i, j) * x[j];
            }
            Eigenvectors.at(i, k) = val;
        }
        
        // Normalize column k
        double norm = 0;
        for(int i=0; i<n; ++i) norm += Eigenvectors.at(i, k) * Eigenvectors.at(i, k);
        norm = std::sqrt(norm);
        for(int i=0; i<n; ++i) Eigenvectors.at(i, k) /= norm;
    }

    return std::make_pair(Eigenvalues, Eigenvectors);
}

std::tuple<Matrix, Matrix, Matrix> Decomposer::SVD(const Matrix& A) {
    // A = U * S * V^T
    // A^T * A = V * S^2 * V^T
    // A * A^T = U * S^2 * U^T
    
    int m = A.getRows();
    int n = A.getCols();
    
    Matrix At = A.transpose();
    Matrix AtA = At * A;
    Logger::getInstance().log("Starting SVD");
    Logger::getInstance().log("Computing Eigenvalues of A^T * A");
    
    auto [Evals, V] = Eigen(AtA);
    
    // Sort eigenvalues and eigenvectors
    std::vector<std::pair<double, int>> evalPairs;
    for (int i = 0; i < n; ++i) {
        evalPairs.push_back({Evals.at(i, i), i});
    }
    std::sort(evalPairs.rbegin(), evalPairs.rend()); // Descending
    
    Matrix S(m, n);
    Matrix V_sorted(n, n);
    Matrix U(m, m);
    
    for (int i = 0; i < n; ++i) {
        double sigma = std::sqrt(std::max(0.0, evalPairs[i].first));
        if (i < m) S.at(i, i) = sigma;
        
        int oldIndex = evalPairs[i].second;
        for (int k = 0; k < n; ++k) {
            V_sorted.at(k, i) = V.at(k, oldIndex);
        }
        
        // Calculate U columns: u_i = A * v_i / sigma
        if (sigma > 1e-10) {
            std::vector<double> v_col(n);
            for(int k=0; k<n; ++k) v_col[k] = V_sorted.at(k, i);
            
            // Multiply A * v_col
            std::vector<double> u_col(m, 0.0);
            for(int r=0; r<m; ++r) {
                for(int c=0; c<n; ++c) {
                    u_col[r] += A.at(r, c) * v_col[c];
                }
            }
            
            for(int k=0; k<m; ++k) U.at(k, i) = u_col[k] / sigma;
        }
    }
    
    // Fill remaining U columns if m > n using Gram-Schmidt or similar if needed
    // For now, this basic implementation assumes we care mostly about the first n singular values/vectors
    
    return std::make_tuple(U, S, V_sorted);
}

std::pair<double, std::vector<double>> Decomposer::PowerIteration(const Matrix& A, int maxIter, double tol) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    
    std::vector<double> b(n, 1.0); // Initial guess
    b = VectorOps::normalize(b);
    double lambda = 0;
    
    Logger::getInstance().log("Starting Power Iteration");
    
    for (int iter = 0; iter < maxIter; ++iter) {
        // b_k+1 = A * b_k
        std::vector<double> b_next(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                b_next[i] += A.at(i, j) * b[j];
            }
        }
        
        double norm = VectorOps::norm(b_next);
        b_next = VectorOps::normalize(b_next);
        
        // Rayleigh quotient: (b_next . A * b_next) / (b_next . b_next) -> approximates lambda
        // Or simpler: norm is approx lambda if A has positive dominant eigenvalue
        // Better: lambda = dot(b_next, A*b_next)
        
        std::vector<double> Ab(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Ab[i] += A.at(i, j) * b_next[j];
            }
        }
        double new_lambda = VectorOps::dot(b_next, Ab);
        
        if (std::abs(new_lambda - lambda) < tol) {
            Logger::getInstance().log("Converged at iteration " + std::to_string(iter));
            return std::make_pair(new_lambda, b_next);
        }
        lambda = new_lambda;
        b = b_next;
    }
    
    return std::make_pair(lambda, b);
}
