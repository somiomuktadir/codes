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
    
    // Iterative refinement for better accuracy
    static std::vector<double> solveRefined(const Matrix& A, const std::vector<double>& b, int maxIter = 5);

    // Least Squares: Solve A * x = b (approx) for overdetermined systems
    static std::vector<double> leastSquares(const Matrix& A, const std::vector<double>& b);
    
    // Solves Ax = b for symmetric positive-definite systems using Cholesky decomposition
    static std::vector<double> solveCholesky(const Matrix& A, const std::vector<double>& b);
    
    // Calculates determinant using Gaussian Elimination to upper triangular form
    static double determinant(const Matrix& A);
    
    // Calculates inverse using Gauss-Jordan Elimination
    static Matrix inverse(const Matrix& A);

    // Matrix power
    static Matrix power(const Matrix& A, int n);

private:
    // Helper to augment matrix A with vector b
    static Matrix augment(const Matrix& A, const std::vector<double>& b);
    
    // Helper to augment matrix A with Identity matrix
    static Matrix augmentIdentity(const Matrix& A);
};

} // namespace LinAlg

#endif // LINEAR_SOLVER_H
