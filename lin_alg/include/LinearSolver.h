#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "Matrix.h"
#include <vector>
#include <cmath>
#include <stdexcept>

namespace LinAlg {

class LinearSolver {
public:
    // Solves Ax = b using Gaussian Elimination with partial pivoting
    static std::vector<double> solve(const Matrix& A, const std::vector<double>& b);
    
    // Solves Ax = b using Cholesky decomposition (requires A to be SPD)
    static std::vector<double> solveCholesky(const Matrix& A, const std::vector<double>& b);
    
    // Least Squares solution for overdetermined systems (Ax = b)
    // Returns x that minimizes ||Ax - b||
    static std::vector<double> leastSquares(const Matrix& A, const std::vector<double>& b);
    
    // Iterative refinement to improve solution accuracy
    static std::vector<double> solveRefined(const Matrix& A, const std::vector<double>& b, int maxIter = 5);
    
    // Generalized Inverse (Moore-Penrose Pseudo-Inverse)
    static Matrix pseudoInverse(const Matrix& A);
    
    // Solve using Pseudo-Inverse (works for any matrix)
    static std::vector<double> solvePseudoInverse(const Matrix& A, const std::vector<double>& b);

    // Determinant (using LU decomposition)
    static double determinant(const Matrix& A);
    
    // Inverse (using LU decomposition)
    static Matrix inverse(const Matrix& A);
    
    // Power of a matrix
    static Matrix power(const Matrix& A, int n);

private:
    // Helper to augment matrix A with vector b
    static Matrix augment(const Matrix& A, const std::vector<double>& b);
    
    // Helper to augment matrix A with Identity matrix
    static Matrix augmentIdentity(const Matrix& A);
};

} // namespace LinAlg

#endif // LINEAR_SOLVER_H
