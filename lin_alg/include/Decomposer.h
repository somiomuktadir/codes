#ifndef DECOMPOSER_H
#define DECOMPOSER_H

#include "Matrix.h"
#include "VectorOps.h"
#include <vector>
#include <cmath>
#include <utility>
#include <tuple>
#include <stdexcept>

namespace LinAlg {

class Decomposer {
public:
    // PLU Decomposition: PA = LU
    // Returns tuple {P, L, U} where P is a permutation matrix
    static std::tuple<Matrix, Matrix, Matrix> PLU(const Matrix& A);

    // LU Decomposition: A = L * U (with partial pivoting)
    // Returns pair {L, U}
    // DEPRECATED: Use PLU instead
    static std::pair<Matrix, Matrix> LU(const Matrix& A);

    // Cholesky Decomposition: A = L * L^T
    // Returns L. Matrix must be symmetric positive-definite.
    static Matrix Cholesky(const Matrix& A);

    // QR Decomposition: A = Q * R
    // Returns pair {Q, R} using Modified Gram-Schmidt
    static std::pair<Matrix, Matrix> QR(const Matrix& A);

    // Eigen Decomposition: A * v = lambda * v
    // Returns pair {Eigenvalues (diagonal matrix), Eigenvectors (columns)}
    // Uses Schur Decomposition
    static std::pair<Matrix, Matrix> Eigen(const Matrix& A, int maxIter = 1000, double tol = 1e-10);

    // Reduction to Hessenberg form: A = Q * H * Q^T
    // Returns pair {Q, H}
    static std::pair<Matrix, Matrix> Hessenberg(const Matrix& A);

    // Schur Decomposition: A = Q * T * Q^T
    // Returns pair {Q, T} where T is upper triangular (or quasi-upper for complex eigenvalues)
    static std::pair<Matrix, Matrix> Schur(const Matrix& A, int maxIter = 1000, double tol = 1e-10);

    // Singular Value Decomposition: A = U * S * V^T
    // Returns tuple {U, S, V}
    // Uses Golub-Kahan Bidiagonalization
    static std::tuple<Matrix, Matrix, Matrix> SVD(const Matrix& A);

    // Bidiagonalization: A = U * B * V^T
    // Returns tuple {U, B, V} where B is bidiagonal
    static std::tuple<Matrix, Matrix, Matrix> Bidiagonalize(const Matrix& A);

    // Golub-Kahan SVD step on bidiagonal matrix
    static std::tuple<Matrix, Matrix, Matrix> GolubKahanSVD(const Matrix& B);
    
    // Power Iteration: Finds dominant eigenvalue and eigenvector
    // Returns pair {eigenvalue, eigenvector}
    static std::pair<double, std::vector<double>> PowerIteration(const Matrix& A, int maxIter = 1000, double tol = 1e-10);
};

} // namespace LinAlg

#endif // DECOMPOSER_H
