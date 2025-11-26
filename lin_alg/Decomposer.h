#ifndef DECOMPOSER_H
#define DECOMPOSER_H

#include "Matrix.h"
#include "VectorOps.h"
#include <vector>
#include <cmath>
#include <utility>
#include <stdexcept>

class Decomposer {
public:
    // LU Decomposition: A = L * U
    // Returns pair {L, U}
    static std::pair<Matrix, Matrix> LU(const Matrix& A);

    // Cholesky Decomposition: A = L * L^T
    // Returns L. Matrix must be symmetric positive-definite.
    static Matrix Cholesky(const Matrix& A);

    // QR Decomposition: A = Q * R
    // Returns pair {Q, R} using Gram-Schmidt
    static std::pair<Matrix, Matrix> QR(const Matrix& A);

    // Eigen Decomposition: A * v = lambda * v
    // Returns pair {Eigenvalues (diagonal matrix), Eigenvectors (columns)}
    // Uses QR Algorithm
    static std::pair<Matrix, Matrix> Eigen(const Matrix& A, int maxIter = 1000, double tol = 1e-10);

    // Singular Value Decomposition: A = U * S * V^T
    // Returns tuple {U, S, V}
    // Uses Eigen decomposition of A^T * A
    static std::tuple<Matrix, Matrix, Matrix> SVD(const Matrix& A);
    
    // Power Iteration: Finds dominant eigenvalue and eigenvector
    // Returns pair {eigenvalue, eigenvector}
    static std::pair<double, std::vector<double>> PowerIteration(const Matrix& A, int maxIter = 1000, double tol = 1e-10);
};

#endif // DECOMPOSER_H
