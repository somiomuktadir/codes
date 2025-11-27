#ifndef ITERATIVE_SOLVERS_H
#define ITERATIVE_SOLVERS_H

#include "Matrix.h"
#include "SparseMatrix.h"
#include <vector>
#include <string>
#include <functional>

namespace LinAlg {

class IterativeSolvers {
public:
    // Conjugate Gradient method for symmetric positive-definite systems
    // Solves Ax = b where A is SPD
    // Returns solution vector x
    static std::vector<double> CG(const Matrix& A, const std::vector<double>& b,
                                   double tol = 1e-6, int maxIter = 1000,
                                   const std::vector<double>& x0 = {});
    
    // Conjugate Gradient for sparse matrices (more efficient)
    static std::vector<double> CG(const SparseMatrixCSR& A, const std::vector<double>& b,
                                   double tol = 1e-6, int maxIter = 1000,
                                   const std::vector<double>& x0 = {});
    
    // Preconditioned Conjugate Gradient with Jacobi preconditioner
    static std::vector<double> PCG(const Matrix& A, const std::vector<double>& b,
                                    double tol = 1e-6, int maxIter = 1000,
                                    const std::vector<double>& x0 = {});
    
    // GMRES (Generalized Minimal Residual) for general systems
    // Solves Ax = b for general (nonsymmetric) matrices
    // restart: restart parameter for GMRES(m)
    static std::vector<double> GMRES(const Matrix& A, const std::vector<double>& b,
                                      int restart = 30, double tol = 1e-6, 
                                      int maxIter = 1000,
                                      const std::vector<double>& x0 = {});
    
    // GMRES for sparse matrices
    static std::vector<double> GMRES(const SparseMatrixCSR& A, const std::vector<double>& b,
                                      int restart = 30, double tol = 1e-6, 
                                      int maxIter = 1000,
                                      const std::vector<double>& x0 = {});
    
private:
    // Helper: Matrix-vector multiplication (dense)
    static std::vector<double> matvec(const Matrix& A, const std::vector<double>& x);
    
    // Helper: Vector operations
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    static double norm(const std::vector<double>& x);
    static std::vector<double> axpy(double alpha, const std::vector<double>& x, 
                                     const std::vector<double>& y); // y = alpha*x + y
    static std::vector<double> scale(double alpha, const std::vector<double>& x);
    static std::vector<double> add(const std::vector<double>& a, const std::vector<double>& b);
    static std::vector<double> subtract(const std::vector<double>& a, const std::vector<double>& b);
    
    // Arnoldi iteration for GMRES
    static void arnoldi(const Matrix& A, std::vector<std::vector<double>>& Q,
                       Matrix& H, int k);
    
    // Arnoldi iteration for sparse GMRES
    static void arnoldi(const SparseMatrixCSR& A, std::vector<std::vector<double>>& Q,
                       Matrix& H, int k);
    
    // Solve least squares problem min ||Hbar*y - beta*e1||
    static std::vector<double> solveLeastSquares(const Matrix& H, const std::vector<double>& rhs);
    
    // Jacobi preconditioner
    static std::vector<double> jacobiPrecondition(const Matrix& A, const std::vector<double>& r);
};

} // namespace LinAlg

#endif // ITERATIVE_SOLVERS_H
