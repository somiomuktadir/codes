#include "IterativeSolvers.h"
#include "VectorOps.h"
#include "Logger.h"
#include "LinearSolver.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace LinAlg {

// ============================================================================
// Helper Functions
// ============================================================================

std::vector<double> IterativeSolvers::matvec(const Matrix& A, const std::vector<double>& x) {
    int n = A.getRows();
    std::vector<double> result(n, 0.0);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < A.getCols(); ++j) {
            result[i] += A.at(i, j) * x[j];
        }
    }
    
    return result;
}

double IterativeSolvers::dot(const std::vector<double>& a, const std::vector<double>& b) {
    return VectorOps::dot(a, b);
}

double IterativeSolvers::norm(const std::vector<double>& x) {
    return VectorOps::norm(x);
}

std::vector<double> IterativeSolvers::axpy(double alpha, const std::vector<double>& x, 
                                            const std::vector<double>& y) {
    std::vector<double> result(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        result[i] = alpha * x[i] + y[i];
    }
    return result;
}

std::vector<double> IterativeSolvers::scale(double alpha, const std::vector<double>& x) {
    std::vector<double> result(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        result[i] = alpha * x[i];
    }
    return result;
}

std::vector<double> IterativeSolvers::add(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<double> IterativeSolvers::subtract(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<double> IterativeSolvers::jacobiPrecondition(const Matrix& A, const std::vector<double>& r) {
    int n = A.getRows();
    std::vector<double> z(n);
    
    for (int i = 0; i < n; ++i) {
        double diag = A.at(i, i);
        if (std::abs(diag) > 1e-10) {
            z[i] = r[i] / diag;
        } else {
            z[i] = r[i];
        }
    }
    
    return z;
}

// ============================================================================
// Conjugate Gradient (CG) - Dense Matrix
// ============================================================================

std::vector<double> IterativeSolvers::CG(const Matrix& A, const std::vector<double>& b,
                                          double tol, int maxIter,
                                          const std::vector<double>& x0) {
    Logger::getInstance().log("[STEP] Starting Conjugate Gradient (CG) solver");
    
    int n = A.getRows();
    
    // Initial guess
    std::vector<double> x = x0.empty() ? std::vector<double>(n, 0.0) : x0;
    
    // r = b - Ax
    std::vector<double> r = subtract(b, matvec(A, x));
    std::vector<double> p = r;
    
    double rsold = dot(r, r);
    double initialResidual = std::sqrt(rsold);
    
    Logger::getInstance().log("[STEP] Initial residual: " + std::to_string(initialResidual));
    
    for (int iter = 0; iter < maxIter; ++iter) {
        std::vector<double> Ap = matvec(A, p);
        double alpha = rsold / dot(p, Ap);
        
        x = axpy(alpha, p, x);       // x = x + alpha*p
        r = axpy(-alpha, Ap, r);     // r = r - alpha*Ap
        
        double rsnew = dot(r, r);
        double residual = std::sqrt(rsnew);
        
        if (iter % 10 == 0 || residual < tol) {
            Logger::getInstance().log("[STEP] Iteration " + std::to_string(iter) + 
                                     ", residual: " + std::to_string(residual));
        }
        
        if (residual < tol) {
            Logger::getInstance().log("[STEP] CG converged in " + std::to_string(iter) + " iterations");
            return x;
        }
        
        double beta = rsnew / rsold;
        p = axpy(beta, p, r);        // p = r + beta*p
        
        rsold = rsnew;
    }
    
    Logger::getInstance().log("[WARNING] CG did not converge within " + std::to_string(maxIter) + " iterations");
    return x;
}

// ============================================================================
// Conjugate Gradient (CG) - Sparse Matrix
// ============================================================================

std::vector<double> IterativeSolvers::CG(const SparseMatrixCSR& A, const std::vector<double>& b,
                                          double tol, int maxIter,
                                          const std::vector<double>& x0) {
    Logger::getInstance().log("[STEP] Starting Conjugate Gradient (CG) solver for sparse matrix");
    
    int n = A.getRows();
    
    // Initial guess
    std::vector<double> x = x0.empty() ? std::vector<double>(n, 0.0) : x0;
    
    // r = b - Ax
    std::vector<double> r = subtract(b, A.multiply(x));
    std::vector<double> p = r;
    
    double rsold = dot(r, r);
    double initialResidual = std::sqrt(rsold);
    
    Logger::getInstance().log("[STEP] Initial residual: " + std::to_string(initialResidual));
    
    for (int iter = 0; iter < maxIter; ++iter) {
        std::vector<double> Ap = A.multiply(p);
        double alpha = rsold / dot(p, Ap);
        
        x = axpy(alpha, p, x);       // x = x + alpha*p
        r = axpy(-alpha, Ap, r);     // r = r - alpha*Ap
        
        double rsnew = dot(r, r);
        double residual = std::sqrt(rsnew);
        
        if (iter % 10 == 0 || residual < tol) {
            Logger::getInstance().log("[STEP] Iteration " + std::to_string(iter) + 
                                     ", residual: " + std::to_string(residual));
        }
        
        if (residual < tol) {
            Logger::getInstance().log("[STEP] CG converged in " + std::to_string(iter) + " iterations");
            return x;
        }
        
        double beta = rsnew / rsold;
        p = axpy(beta, p, r);        // p = r + beta*p
        
        rsold = rsnew;
    }
    
    Logger::getInstance().log("[WARNING] CG did not converge within " + std::to_string(maxIter) + " iterations");
    return x;
}

// ============================================================================
// Preconditioned Conjugate Gradient (PCG)
// ============================================================================

std::vector<double> IterativeSolvers::PCG(const Matrix& A, const std::vector<double>& b,
                                           double tol, int maxIter,
                                           const std::vector<double>& x0) {
    Logger::getInstance().log("[STEP] Starting Preconditioned Conjugate Gradient (PCG) solver");
    
    int n = A.getRows();
    
    // Initial guess
    std::vector<double> x = x0.empty() ? std::vector<double>(n, 0.0) : x0;
    
    // r = b - Ax
    std::vector<double> r = subtract(b, matvec(A, x));
    std::vector<double> z = jacobiPrecondition(A, r);
    std::vector<double> p = z;
    
    double rzold = dot(r, z);
    double initialResidual = norm(r);
    
    Logger::getInstance().log("[STEP] Initial residual: " + std::to_string(initialResidual));
    
    for (int iter = 0; iter < maxIter; ++iter) {
        std::vector<double> Ap = matvec(A, p);
        double alpha = rzold / dot(p, Ap);
        
        x = axpy(alpha, p, x);       // x = x + alpha*p
        r = axpy(-alpha, Ap, r);     // r = r - alpha*Ap
        
        double residual = norm(r);
        
        if (iter % 10 == 0 || residual < tol) {
            Logger::getInstance().log("[STEP] Iteration " + std::to_string(iter) + 
                                     ", residual: " + std::to_string(residual));
        }
        
        if (residual < tol) {
            Logger::getInstance().log("[STEP] PCG converged in " + std::to_string(iter) + " iterations");
            return x;
        }
        
        z = jacobiPrecondition(A, r);
        double rznew = dot(r, z);
        double beta = rznew / rzold;
        p = axpy(beta, p, z);        // p = z + beta*p
        
        rzold = rznew;
    }
    
    Logger::getInstance().log("[WARNING] PCG did not converge within " + std::to_string(maxIter) + " iterations");
    return x;
}

// ============================================================================
// Arnoldi Iteration
// ============================================================================

void IterativeSolvers::arnoldi(const Matrix& A, std::vector<std::vector<double>>& Q,
                                Matrix& H, int k) {
    int n = A.getRows();
    
    // q = A * Q[k]
    std::vector<double> q = matvec(A, Q[k]);
    
    // Modified Gram-Schmidt orthogonalization
    for (int j = 0; j <= k; ++j) {
        double h_jk = dot(Q[j], q);
        H.at(j, k) = h_jk;
        q = axpy(-h_jk, Q[j], q);
    }
    
    double h_kp1_k = norm(q);
    H.at(k + 1, k) = h_kp1_k;
    
    if (h_kp1_k > 1e-10) {
        Q.push_back(scale(1.0 / h_kp1_k, q));
    } else {
        Q.push_back(std::vector<double>(n, 0.0));
    }
}

void IterativeSolvers::arnoldi(const SparseMatrixCSR& A, std::vector<std::vector<double>>& Q,
                                Matrix& H, int k) {
    int n = A.getRows();
    
    // q = A * Q[k]
    std::vector<double> q = A.multiply(Q[k]);
    
    // Modified Gram-Schmidt orthogonalization
    for (int j = 0; j <= k; ++j) {
        double h_jk = dot(Q[j], q);
        H.at(j, k) = h_jk;
        q = axpy(-h_jk, Q[j], q);
    }
    
    double h_kp1_k = norm(q);
    H.at(k + 1, k) = h_kp1_k;
    
    if (h_kp1_k > 1e-10) {
        Q.push_back(scale(1.0 / h_kp1_k, q));
    } else {
        Q.push_back(std::vector<double>(n, 0.0));
    }
}

// ============================================================================
// GMRES Helper
// ============================================================================

std::vector<double> IterativeSolvers::solveLeastSquares(const Matrix& H, const std::vector<double>& rhs) {
    // Use existing QR solver for least squares
    Matrix A(rhs.size(), H.getCols());
    for (int i = 0; i < (int)rhs.size(); ++i) {
        for (int j = 0; j < H.getCols(); ++j) {
            A.at(i, j) = H.at(i, j);
        }
    }
    
    return LinearSolver::leastSquares(A, rhs);
}

// ============================================================================
// GMRES - Dense Matrix
// ============================================================================

std::vector<double> IterativeSolvers::GMRES(const Matrix& A, const std::vector<double>& b,
                                             int restart, double tol, int maxIter,
                                             const std::vector<double>& x0) {
    Logger::getInstance().log("[STEP] Starting GMRES(" + std::to_string(restart) + ") solver");
    
    int n = A.getRows();
    std::vector<double> x = x0.empty() ? std::vector<double>(n, 0.0) : x0;
    
    double initialResidual = norm(subtract(b, matvec(A, x)));
    Logger::getInstance().log("[STEP] Initial residual: " + std::to_string(initialResidual));
    
    int totalIter = 0;
    
    while (totalIter < maxIter) {
        // r = b - Ax
        std::vector<double> r = subtract(b, matvec(A, x));
        double beta = norm(r);
        
        if (beta < tol) {
            Logger::getInstance().log("[STEP] GMRES converged, residual: " + std::to_string(beta));
            return x;
        }
        
        // Initialize Krylov subspace
        std::vector<std::vector<double>> Q;
        Q.push_back(scale(1.0 / beta, r));
        
        // Upper Hessenberg matrix
        Matrix H(restart + 1, restart, 0.0);
        
        // RHS for least squares
        std::vector<double> g(restart + 1, 0.0);
        g[0] = beta;
        
        int m = std::min(restart, maxIter - totalIter);
        
        for (int k = 0; k < m; ++k) {
            totalIter++;
            arnoldi(A, Q, H, k);
            
            // Solve least squares problem
            std::vector<double> y = solveLeastSquares(H, g);
            
            // Compute residual norm
            std::vector<double> Hy(k + 2, 0.0);
            for (int i = 0; i <= k + 1; ++i) {
                for (int j = 0; j <= k; ++j) {
                    Hy[i] += H.at(i, j) * y[j];
                }
            }
            double residual = norm(subtract(g, Hy));
            
            if (k % 5 == 0 || residual < tol) {
                Logger::getInstance().log("[STEP] GMRES iteration " + std::to_string(totalIter) + 
                                         ", residual: " + std::to_string(residual));
            }
            
            if (residual < tol) {
                // Compute solution
                for (int j = 0; j <= k; ++j) {
                    x = axpy(y[j], Q[j], x);
                }
                Logger::getInstance().log("[STEP] GMRES converged in " + std::to_string(totalIter) + " iterations");
                return x;
            }
        }
        
        // Update solution with current Krylov subspace
        std::vector<double> y = solveLeastSquares(H, g);
        for (int j = 0; j < m; ++j) {
            x = axpy(y[j], Q[j], x);
        }
    }
    
    Logger::getInstance().log("[WARNING] GMRES did not converge within " + std::to_string(maxIter) + " iterations");
    return x;
}

// ============================================================================
// GMRES - Sparse Matrix
// ============================================================================

std::vector<double> IterativeSolvers::GMRES(const SparseMatrixCSR& A, const std::vector<double>& b,
                                             int restart, double tol, int maxIter,
                                             const std::vector<double>& x0) {
    Logger::getInstance().log("[STEP] Starting GMRES(" + std::to_string(restart) + ") solver for sparse matrix");
    
    int n = A.getRows();
    std::vector<double> x = x0.empty() ? std::vector<double>(n, 0.0) : x0;
    
    double initialResidual = norm(subtract(b, A.multiply(x)));
    Logger::getInstance().log("[STEP] Initial residual: " + std::to_string(initialResidual));
    
    int totalIter = 0;
    
    while (totalIter < maxIter) {
        // r = b - Ax
        std::vector<double> r = subtract(b, A.multiply(x));
        double beta = norm(r);
        
        if (beta < tol) {
            Logger::getInstance().log("[STEP] GMRES converged, residual: " + std::to_string(beta));
            return x;
        }
        
        // Initialize Krylov subspace
        std::vector<std::vector<double>> Q;
        Q.push_back(scale(1.0 / beta, r));
        
        // Upper Hessenberg matrix
        Matrix H(restart + 1, restart, 0.0);
        
        // RHS for least squares
        std::vector<double> g(restart + 1, 0.0);
        g[0] = beta;
        
        int m = std::min(restart, maxIter - totalIter);
        
        for (int k = 0; k < m; ++k) {
            totalIter++;
            arnoldi(A, Q, H, k);
            
            // Solve least squares problem
            std::vector<double> y = solveLeastSquares(H, g);
            
            // Compute residual norm
            std::vector<double> Hy(k + 2, 0.0);
            for (int i = 0; i <= k + 1; ++i) {
                for (int j = 0; j <= k; ++j) {
                    Hy[i] += H.at(i, j) * y[j];
                }
            }
            double residual = norm(subtract(g, Hy));
            
            if (k % 5 == 0 || residual < tol) {
                Logger::getInstance().log("[STEP] GMRES iteration " + std::to_string(totalIter) + 
                                         ", residual: " + std::to_string(residual));
            }
            
            if (residual < tol) {
                // Compute solution
                for (int j = 0; j <= k; ++j) {
                    x = axpy(y[j], Q[j], x);
                }
                Logger::getInstance().log("[STEP] GMRES converged in " + std::to_string(totalIter) + " iterations");
                return x;
            }
        }
        
        // Update solution with current Krylov subspace
        std::vector<double> y = solveLeastSquares(H, g);
        for (int j = 0; j < m; ++j) {
            x = axpy(y[j], Q[j], x);
        }
    }
    
    Logger::getInstance().log("[WARNING] GMRES did not converge within " + std::to_string(maxIter) + " iterations");
    return x;
}

} // namespace LinAlg
