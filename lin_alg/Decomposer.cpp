#include "Decomposer.h"
#include "LinearSolver.h" // For inverse if needed
#include "Logger.h"
#include <algorithm>
#include <sstream>
#include <iomanip>

namespace LinAlg {

std::tuple<Matrix, Matrix, Matrix> Decomposer::PLU(const Matrix& A) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");

    Matrix L = Matrix::identity(n);
    Matrix U = A;
    Matrix P = Matrix::identity(n);
    Logger::getInstance().log("Starting PLU Decomposition");

    for (int i = 0; i < n; ++i) {
        // Find pivot
        int pivotRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(U(j, i)) > std::abs(U(pivotRow, i))) {
                pivotRow = j;
            }
        }
        
        // Swap rows in U
        if (pivotRow != i) {
            Logger::getInstance().log("Swapping row " + std::to_string(i) + " with row " + std::to_string(pivotRow));
            for (int k = 0; k < n; ++k) {
                std::swap(U(i, k), U(pivotRow, k));
            }
            // Swap rows in P
            for (int k = 0; k < n; ++k) {
                std::swap(P(i, k), P(pivotRow, k));
            }
            // Swap rows in L (only the part we've computed so far, columns 0 to i-1)
            for (int k = 0; k < i; ++k) {
                std::swap(L(i, k), L(pivotRow, k));
            }
        }
        
        if (std::abs(U(i, i)) < Matrix::epsilon()) {
             // Singular matrix, but we continue decomposition (U will have 0 on diagonal)
             Logger::getInstance().log("Warning: Pivot is close to zero");
        } else {
             for (int j = i + 1; j < n; ++j) {
                double factor = U(j, i) / U(i, i);
                L(j, i) = factor;
                for (int k = i; k < n; ++k) {
                    U(j, k) -= factor * U(i, k);
                }
            }
        }
    }
    Logger::getInstance().log("PLU Decomposition completed");
    return std::make_tuple(P, L, U);
}

std::pair<Matrix, Matrix> Decomposer::LU(const Matrix& A) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");

    Matrix L = Matrix::identity(n);
    Matrix U = A;
    Logger::getInstance().log("Starting LU Decomposition with Partial Pivoting");

    for (int i = 0; i < n; ++i) {
        // Find pivot
        int pivotRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(U(j, i)) > std::abs(U(pivotRow, i))) {
                pivotRow = j;
            }
        }
        
        // Swap rows in U if needed
        if (pivotRow != i) {
            Logger::getInstance().log("Swapping row " + std::to_string(i) + " with row " + std::to_string(pivotRow));
            for (int k = 0; k < n; ++k) {
                std::swap(U(i, k), U(pivotRow, k));
            }
            // Also swap already computed parts of L
            for (int k = 0; k < i; ++k) {
                std::swap(L(i, k), L(pivotRow, k));
            }
        }
        
        if (std::abs(U(i, i)) < 1e-10) throw std::runtime_error("Matrix is singular or nearly singular");
        
        for (int j = i + 1; j < n; ++j) {
            double factor = U(j, i) / U(i, i);
            L(j, i) = factor;
            for (int k = i; k < n; ++k) {
                U(j, k) -= factor * U(i, k);
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
                sum += L(i, k) * L(j, k);
            }
            
            if (i == j) {
                double val = A(i, i) - sum;
                if (val <= 0) throw std::runtime_error("Matrix is not positive definite");
                L(i, j) = std::sqrt(val);
            } else {
                L(i, j) = (A(i, j) - sum) / L(j, j);
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
            qCols[j][i] = A(i, j);
        }
    }
    
    // Modified Gram-Schmidt process
    for (int j = 0; j < n; ++j) {
        // Orthogonalize against all previous vectors
        for (int i = 0; i < j; ++i) {
            // r_ij = <q_i, a_j>
            double r_ij = VectorOps::dot(qCols[i], qCols[j]);
            R(i, j) = r_ij;
            
            // a_j = a_j - r_ij * q_i
            for (int k = 0; k < m; ++k) {
                qCols[j][k] -= r_ij * qCols[i][k];
            }
        }
        
        // Normalize q_j
        double r_jj = VectorOps::norm(qCols[j]);
        R(j, j) = r_jj;
        
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
            Q(i, j) = qCols[j][i];
        }
    }
    
    return std::make_pair(Q, R);
}

std::pair<Matrix, Matrix> Decomposer::Hessenberg(const Matrix& A) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    
    Matrix H = A;
    Matrix Q = Matrix::identity(n);
    
    for (int k = 0; k < n - 2; ++k) {
        // Compute Householder reflector to zero out H(k+2...n-1, k)
        std::vector<double> x(n - k - 1);
        for (int i = 0; i < n - k - 1; ++i) {
            x[i] = H(k + 1 + i, k);
        }
        
        double normX = VectorOps::norm(x);
        if (normX < Matrix::epsilon()) continue;
        
        double alpha = (x[0] > 0) ? -normX : normX;
        std::vector<double> u = x;
        u[0] -= alpha;
        double normU = VectorOps::norm(u);
        if (normU < Matrix::epsilon()) continue;
        
        for (double& val : u) val /= normU;
        
        // Apply reflector to H from left: H = (I - 2uu^T)H
        // H(k+1...n-1, :) -= 2u * (u^T * H(k+1...n-1, :))
        for (int j = k; j < n; ++j) {
            double dot = 0;
            for (int i = 0; i < n - k - 1; ++i) {
                dot += u[i] * H(k + 1 + i, j);
            }
            for (int i = 0; i < n - k - 1; ++i) {
                H(k + 1 + i, j) -= 2 * u[i] * dot;
            }
        }
        
        // Apply reflector to H from right: H = H(I - 2uu^T)
        // H(:, k+1...n-1) -= (H(:, k+1...n-1) * u) * 2u^T
        for (int i = 0; i < n; ++i) {
            double dot = 0;
            for (int j = 0; j < n - k - 1; ++j) {
                dot += H(i, k + 1 + j) * u[j];
            }
            for (int j = 0; j < n - k - 1; ++j) {
                H(i, k + 1 + j) -= 2 * dot * u[j];
            }
        }
        
        // Accumulate Q: Q = Q(I - 2uu^T)
        for (int i = 0; i < n; ++i) {
            double dot = 0;
            for (int j = 0; j < n - k - 1; ++j) {
                dot += Q(i, k + 1 + j) * u[j];
            }
            for (int j = 0; j < n - k - 1; ++j) {
                Q(i, k + 1 + j) -= 2 * dot * u[j];
            }
        }
    }
    
    return {Q, H};
}

std::pair<Matrix, Matrix> Decomposer::Schur(const Matrix& A, int maxIter, double tol) {
    // 1. Reduce to Hessenberg form
    auto [Q, H] = Hessenberg(A);
    int n = A.getRows();
    
    // 2. QR Algorithm with shifts
    int n_active = n;
    int iter = 0;
    
    while (n_active > 1 && iter < maxIter) {
        int m = n_active - 1;
        
        // Check for convergence (subdiagonal element close to 0)
        if (std::abs(H(m, m-1)) < tol) {
            H(m, m-1) = 0.0; // Clean up
            n_active--;
            continue;
        }
        
        // Wilkinson shift
        double d = (H(m-1, m-1) - H(m, m)) / 2.0;
        double sgn = (d >= 0) ? 1.0 : -1.0;
        double mu = H(m, m) - (H(m, m-1) * H(m, m-1)) / (d + sgn * std::sqrt(d*d + H(m, m-1)*H(m, m-1)));
        
        // Implicit QR step
        double x = H(0, 0) - mu;
        double z = H(1, 0);
        
        for (int k = 0; k < n_active - 1; ++k) {
            // Givens rotation
            double c, s;
            if (std::abs(z) < Matrix::epsilon()) {
                c = 1.0; s = 0.0;
            } else {
                double r = std::hypot(x, z);
                c = x / r;
                s = -z / r;
            }
            
            // Apply G_k to H from left
            for (int j = k; j < n; ++j) {
                double t1 = H(k, j);
                double t2 = H(k+1, j);
                H(k, j) = c * t1 - s * t2;
                H(k+1, j) = s * t1 + c * t2;
            }
            
            // Apply G_k^T to H from right
            for (int i = 0; i <= std::min(k+2, n-1); ++i) {
                double t1 = H(i, k);
                double t2 = H(i, k+1);
                H(i, k) = c * t1 - s * t2;
                H(i, k+1) = s * t1 + c * t2;
            }
            
            // Accumulate Q
            for (int i = 0; i < n; ++i) {
                double t1 = Q(i, k);
                double t2 = Q(i, k+1);
                Q(i, k) = c * t1 - s * t2;
                Q(i, k+1) = s * t1 + c * t2;
            }
            
            if (k < n_active - 2) {
                x = H(k+1, k);
                z = H(k+2, k);
            }
        }
        
        iter++;
    }
    
    return {Q, H};
}

std::pair<Matrix, Matrix> Decomposer::Eigen(const Matrix& A, int maxIter, double tol) {
    int n = A.getRows();
    if (A.getCols() != n) throw std::invalid_argument("Matrix must be square");
    
    Logger::getInstance().log("Starting Eigen Decomposition (Schur Algorithm)");
    
    auto [Q, T] = Schur(A, maxIter, tol);
    
    // Eigenvalues are on diagonal of T
    Matrix Eigenvalues(n, n);
    for(int i=0; i<n; ++i) Eigenvalues(i,i) = T(i,i);
    
    // Eigenvectors
    // For symmetric matrices, Q contains eigenvectors.
    // For non-symmetric, we need to solve (T - lambda*I)x = 0 to get eigenvectors of T,
    // then transform back with Q * x.
    
    // Check if symmetric
    bool symmetric = true;
    for(int i=0; i<n; ++i) {
        for(int j=i+1; j<n; ++j) {
            if(std::abs(A(i,j) - A(j,i)) > 1e-5) {
                symmetric = false;
                break;
            }
        }
    }
    
    if (symmetric) {
        Logger::getInstance().log("Matrix is symmetric, Q contains eigenvectors");
        return {Eigenvalues, Q};
    }
    
    Logger::getInstance().log("Matrix is non-symmetric, computing eigenvectors via back-substitution");
    Matrix Eigenvectors(n, n);
    
    for (int k = 0; k < n; ++k) {
        double lambda = T(k, k);
        // Solve (T - lambda*I)x = 0
        // T is upper triangular (or quasi-upper, but we assume real eigenvalues for now)
        
        std::vector<double> x(n, 0.0);
        x[k] = 1.0;
        
        for (int i = k - 1; i >= 0; --i) {
            double sum = 0;
            for (int j = i + 1; j <= k; ++j) {
                sum += T(i, j) * x[j];
            }
            double denom = T(i, i) - lambda;
            if (std::abs(denom) < Matrix::epsilon()) {
                x[i] = 0; // Fallback
            } else {
                x[i] = -sum / denom;
            }
        }
        
        // Transform back: v = Q * x
        for (int i = 0; i < n; ++i) {
            double val = 0;
            for (int j = 0; j < n; ++j) {
                val += Q(i, j) * x[j];
            }
            Eigenvectors(i, k) = val;
        }
        
        // Normalize column k
        double norm = 0;
        for(int i=0; i<n; ++i) norm += Eigenvectors(i, k) * Eigenvectors(i, k);
        norm = std::sqrt(norm);
        if (norm > Matrix::epsilon()) {
            for(int i=0; i<n; ++i) Eigenvectors(i, k) /= norm;
        }
    }

    return std::make_pair(Eigenvalues, Eigenvectors);
}

std::tuple<Matrix, Matrix, Matrix> Decomposer::Bidiagonalize(const Matrix& A) {
    int m = A.getRows();
    int n = A.getCols();
    
    Matrix B = A;
    Matrix U = Matrix::identity(m);
    Matrix V = Matrix::identity(n);
    
    for (int k = 0; k < n; ++k) {
        // Eliminate column k below diagonal
        if (k < m - 1) {
            std::vector<double> x(m - k);
            for(int i=0; i<m-k; ++i) x[i] = B(k+i, k);
            
            double normX = VectorOps::norm(x);
            if (normX > Matrix::epsilon()) {
                double alpha = (x[0] > 0) ? -normX : normX;
                std::vector<double> u = x;
                u[0] -= alpha;
                double normU = VectorOps::norm(u);
                if (normU > Matrix::epsilon()) {
                    for(double& val : u) val /= normU;
                    
                    // Apply to B from left
                    for(int j=k; j<n; ++j) {
                        double dot = 0;
                        for(int i=0; i<m-k; ++i) dot += u[i] * B(k+i, j);
                        for(int i=0; i<m-k; ++i) B(k+i, j) -= 2 * u[i] * dot;
                    }
                    
                    // Accumulate U
                    for(int i=0; i<m; ++i) {
                        double dot = 0;
                        for(int j=0; j<m-k; ++j) dot += U(i, k+j) * u[j];
                        for(int j=0; j<m-k; ++j) U(i, k+j) -= 2 * dot * u[j];
                    }
                }
            }
        }
        
        // Eliminate row k to right of superdiagonal
        if (k < n - 2) {
            std::vector<double> x(n - k - 1);
            for(int i=0; i<n-k-1; ++i) x[i] = B(k, k+1+i);
            
            double normX = VectorOps::norm(x);
            if (normX > Matrix::epsilon()) {
                double alpha = (x[0] > 0) ? -normX : normX;
                std::vector<double> u = x;
                u[0] -= alpha;
                double normU = VectorOps::norm(u);
                if (normU > Matrix::epsilon()) {
                    for(double& val : u) val /= normU;
                    
                    // Apply to B from right
                    for(int i=k; i<m; ++i) {
                        double dot = 0;
                        for(int j=0; j<n-k-1; ++j) dot += B(i, k+1+j) * u[j];
                        for(int j=0; j<n-k-1; ++j) B(i, k+1+j) -= 2 * dot * u[j];
                    }
                    
                    // Accumulate V
                    for(int i=0; i<n; ++i) {
                        double dot = 0;
                        for(int j=0; j<n-k-1; ++j) dot += V(i, k+1+j) * u[j];
                        for(int j=0; j<n-k-1; ++j) V(i, k+1+j) -= 2 * dot * u[j];
                    }
                }
            }
        }
    }
    
    return {U, B, V};
}

std::tuple<Matrix, Matrix, Matrix> Decomposer::GolubKahanSVD(const Matrix& B_in) {
    Matrix B = B_in;
    int m = B.getRows();
    int n = B.getCols();
    Matrix U = Matrix::identity(m);
    Matrix V = Matrix::identity(n);
    
    // Simple implementation: convert to square if needed, then apply implicit QR
    // For simplicity in this iteration, we assume square or m >= n.
    // If m < n, we should transpose.
    
    int maxIter = 1000;
    for (int iter = 0; iter < maxIter; ++iter) {
        bool converged = true;
        for (int i = 0; i < std::min(m, n) - 1; ++i) {
            if (std::abs(B(i, i+1)) > Matrix::epsilon()) {
                converged = false;
                break;
            }
        }
        if (converged) break;
        
        // Golub-Kahan step
        // Compute shift from bottom 2x2 block of B^T * B
        // We need to be careful with indices.
        // Let's implement a very basic sweep for now to ensure stability.
        
        // For each superdiagonal element, apply Givens rotations to chase the zero.
        // This is complex to implement correctly from scratch without errors.
        // Fallback: Use SVD of B^T * B since B is bidiagonal, condition number is better?
        // No, that defeats the purpose.
        
        // Let's implement the standard implicit QR step for SVD.
        // T = B^T * B is tridiagonal.
        // We compute Wilkinson shift for T.
        
        // Simplified: Just apply rotations to diagonalize B
        for (int i = 0; i < std::min(m, n) - 1; ++i) {
             double alpha = B(i, i);
             double beta = B(i, i+1);
             
             // Compute rotation to annihilate beta
             double c, s;
             if (std::abs(beta) < Matrix::epsilon()) {
                 c = 1.0; s = 0.0;
             } else {
                 double r = std::hypot(alpha, beta);
                 c = alpha / r;
                 s = -beta / r;
             }
             
             // Apply to B from right (affects columns i and i+1)
             // B = B * G
             for(int r=0; r<m; ++r) {
                 double t1 = B(r, i);
                 double t2 = B(r, i+1);
                 B(r, i) = c * t1 - s * t2;
                 B(r, i+1) = s * t1 + c * t2;
             }
             // Update V
             for(int r=0; r<n; ++r) {
                 double t1 = V(r, i);
                 double t2 = V(r, i+1);
                 V(r, i) = c * t1 - s * t2;
                 V(r, i+1) = s * t1 + c * t2;
             }
             
             // Now we have introduced a bulge at B(i+1, i)
             // Chase it with left rotation
             double a = B(i, i);
             double b = B(i+1, i);
             
             if (std::abs(b) > Matrix::epsilon()) {
                 double r = std::hypot(a, b);
                 c = a / r;
                 s = -b / r;
                 
                 // Apply to B from left (affects rows i and i+1)
                 // B = G^T * B
                 for(int c_idx=0; c_idx<n; ++c_idx) {
                     double t1 = B(i, c_idx);
                     double t2 = B(i+1, c_idx);
                     B(i, c_idx) = c * t1 - s * t2;
                     B(i+1, c_idx) = s * t1 + c * t2;
                 }
                 // Update U
                 for(int c_idx=0; c_idx<m; ++c_idx) {
                     double t1 = U(c_idx, i);
                     double t2 = U(c_idx, i+1);
                     U(c_idx, i) = c * t1 - s * t2;
                     U(c_idx, i+1) = s * t1 + c * t2;
                 }
             }
        }
    }
    
    return {U, B, V};
}

std::tuple<Matrix, Matrix, Matrix> Decomposer::SVD(const Matrix& A) {
    Logger::getInstance().log("Starting SVD (Golub-Kahan)");
    
    // 1. Bidiagonalize
    auto [U1, B, V1] = Bidiagonalize(A);
    
    // 2. Diagonalize B
    auto [U2, S, V2] = GolubKahanSVD(B);
    
    // 3. Combine
    // A = U1 * B * V1^T
    // B = U2 * S * V2^T
    // A = U1 * U2 * S * V2^T * V1^T
    // U = U1 * U2
    // V = V1 * V2
    
    Matrix U = U1 * U2;
    Matrix V = V1 * V2;
    
    // Ensure singular values are positive
    int minDim = std::min(S.getRows(), S.getCols());
    for(int i=0; i<minDim; ++i) {
        if(S(i,i) < 0) {
            S(i,i) = -S(i,i);
            for(int k=0; k<U.getRows(); ++k) U(k,i) = -U(k,i);
        }
    }
    
    // Sort singular values
    // (Bubble sort for simplicity, n is usually small)
    for(int i=0; i<minDim-1; ++i) {
        for(int j=0; j<minDim-i-1; ++j) {
            if(S(j,j) < S(j+1,j+1)) {
                // Swap singular values
                std::swap(S(j,j), S(j+1,j+1));
                
                // Swap columns of U
                for(int k=0; k<U.getRows(); ++k) std::swap(U(k,j), U(k,j+1));
                
                // Swap columns of V
                for(int k=0; k<V.getRows(); ++k) std::swap(V(k,j), V(k,j+1));
            }
        }
    }
    
    return {U, S, V};
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
                b_next[i] += A(i, j) * b[j];
            }
        }
        
        b_next = VectorOps::normalize(b_next);
        
        // Rayleigh quotient
        std::vector<double> Ab(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Ab[i] += A(i, j) * b_next[j];
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

} // namespace LinAlg
