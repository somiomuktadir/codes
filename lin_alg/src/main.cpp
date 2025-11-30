#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include <complex>
#include "Matrix.h"
#include "VectorOps.h"
#include "LinearSolver.h"
#include "Decomposer.h"
#include "Analysis.h"
#include "Logger.h"
#include "SparseMatrix.h"
#include "IterativeSolvers.h"
#include "Statistics.h"
#include "ComplexMatrix.h"
#include "QuadraticForm.h"

using namespace LinAlg;

void clearInput() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

Matrix inputMatrix(const std::string& name) {
    int r, c;
    std::cout << "Enter dimensions for " << name << " (rows cols): ";
    std::cin >> r >> c;
    Matrix M(r, c);
    std::cout << "Enter elements row by row:" << std::endl;
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            std::cin >> M(i, j);
        }
    }
    return M;
}

std::vector<double> inputVector(const std::string& name) {
    int n;
    std::cout << "Enter size for " << name << ": ";
    std::cin >> n;
    std::vector<double> v(n);
    std::cout << "Enter elements: ";
    for (int i = 0; i < n; ++i) {
        std::cin >> v[i];
    }
    return v;
}

void printMenu() {
    std::cout << "\n=== Linear Algebra CLI ===";
    #if defined(USE_BLAS) || defined(_OPENMP)
    std::cout << " [";
    #ifdef USE_BLAS
    std::cout << "BLAS";
    #endif
    #ifdef _OPENMP
    #ifdef USE_BLAS
    std::cout << "+";
    #endif
    std::cout << "OpenMP";
    #endif
    std::cout << "]";
    #endif
    std::cout << std::endl;
    
    std::cout << "1. Matrix Operations (+, -, *, T, Norms, Kronecker)" << std::endl;
    std::cout << "2. Vector Operations (Dot, Norm, Cross)" << std::endl;
    std::cout << "3. Linear Solver (Ax = b)" << std::endl;
    std::cout << "4. Determinant, Inverse & Pseudo-Inverse" << std::endl;
    std::cout << "5. Decompositions (LU, Cholesky, QR, Eigen, SVD, Diagonalization)" << std::endl;
    std::cout << "6. Analysis (PCA, Rank, Trace, Quadratic Form)" << std::endl;
    std::cout << "7. Statistics (Mean, Variance, Covariance, Correlation)" << std::endl;
    std::cout << "8. Matrix Manipulation (Submatrix, Stack, Hadamard)" << std::endl;
    std::cout << "9. Toggle Verbose Mode (Current: " << (Logger::getInstance().isEnabled() ? "ON" : "OFF") << ")" << std::endl;
    std::cout << "10. Sparse Matrices (CSR/CSC Format)" << std::endl;
    std::cout << "11. Iterative Solvers (CG, GMRES)" << std::endl;
    std::cout << "12. Complex Matrices" << std::endl;
    std::cout << "0. Exit" << std::endl;
    std::cout << "Select option: ";
}

int main() {
    int choice;
    while (true) {
        printMenu();
        if (!(std::cin >> choice)) {
            clearInput();
            continue;
        }

        try {
            if (choice == 0) break;
            
            if (choice == 9) {
                if (Logger::getInstance().isEnabled()) Logger::getInstance().disable();
                else Logger::getInstance().enable();
                continue;
            }
            
            switch (choice) {
                case 1: { // Matrix Ops
                    std::cout << "1. Add  2. Subtract  3. Multiply  4. Transpose  5. Power  6. Hadamard  7. Norms  8. Kronecker" << std::endl;
                    int sub; std::cin >> sub;
                    if (sub == 4) {
                        Matrix A = inputMatrix("Matrix A");
                        std::cout << "Transpose:" << std::endl;
                        A.transpose().print();
                    } else if (sub == 5) {
                        Matrix A = inputMatrix("Matrix A");
                        int n;
                        std::cout << "Enter power n: "; std::cin >> n;
                        std::cout << "A^" << n << ":" << std::endl;
                        LinearSolver::power(A, n).print();
                    } else if (sub == 6) {
                        Matrix A = inputMatrix("Matrix A");
                        Matrix B = inputMatrix("Matrix B");
                        std::cout << "Hadamard product (element-wise):" << std::endl;
                        A.hadamard(B).print();
                    } else if (sub == 7) {
                        Matrix A = inputMatrix("Matrix A");
                        std::cout << "Frobenius Norm: " << A.frobeniusNorm() << std::endl;
                        std::cout << "L1 Norm: " << A.l1Norm() << std::endl;
                        std::cout << "L-Inf Norm: " << A.lInfNorm() << std::endl;
                    } else if (sub == 8) {
                        Matrix A = inputMatrix("Matrix A");
                        Matrix B = inputMatrix("Matrix B");
                        std::cout << "Kronecker Product:" << std::endl;
                        Matrix::kroneckerProduct(A, B).print();
                    } else {
                        Matrix A = inputMatrix("Matrix A");
                        Matrix B = inputMatrix("Matrix B");
                        if (sub == 1) (A + B).print();
                        else if (sub == 2) (A - B).print();
                        else if (sub == 3) (A * B).print();
                    }
                    break;
                }
                case 2: { // Vector Ops
                    std::cout << "1. Dot  2. Norm  3. Cross" << std::endl;
                    int sub; std::cin >> sub;
                    if (sub == 2) {
                        std::vector<double> v = inputVector("Vector");
                        std::cout << "Norm: " << VectorOps::norm(v) << std::endl;
                    } else {
                        std::vector<double> v1 = inputVector("Vector 1");
                        std::vector<double> v2 = inputVector("Vector 2");
                        if (sub == 1) std::cout << "Dot: " << VectorOps::dot(v1, v2) << std::endl;
                        else if (sub == 3) {
                            std::cout << "Cross: ";
                            VectorOps::print(VectorOps::cross(v1, v2));
                        }
                    }
                    break;
                }
                case 3: { // Linear Solver
                    std::cout << "1. Standard  2. Cholesky (for SPD)  3. Least Squares  4. Refined" << std::endl;
                    int sub; std::cin >> sub;
                    Matrix A = inputMatrix("Matrix A");
                    std::vector<double> b = inputVector("Vector b");
                    std::vector<double> x;
                    
                    if (sub == 1) x = LinearSolver::solve(A, b);
                    else if (sub == 2) x = LinearSolver::solveCholesky(A, b);
                    else if (sub == 3) x = LinearSolver::leastSquares(A, b);
                    else if (sub == 4) x = LinearSolver::solveRefined(A, b);
                    
                    std::cout << "Solution x: ";
                    VectorOps::print(x);
                    break;
                }
                case 4: { // Det & Inv
                    std::cout << "1. Determinant  2. Inverse  3. Pseudo-Inverse" << std::endl;
                    int sub; std::cin >> sub;
                    Matrix A = inputMatrix("Matrix A");
                    
                    if (sub == 1) {
                        std::cout << "Determinant: " << LinearSolver::determinant(A) << std::endl;
                    } else if (sub == 2) {
                        try {
                            std::cout << "Inverse:" << std::endl;
                            LinearSolver::inverse(A).print();
                        } catch (const std::exception& e) {
                            std::cout << "Inverse not possible: " << e.what() << std::endl;
                        }
                    } else if (sub == 3) {
                        std::cout << "Pseudo-Inverse:" << std::endl;
                        LinearSolver::pseudoInverse(A).print();
                    }
                    break;
                }
                case 5: { // Decompositions
                    std::cout << "1. LU  2. Cholesky  3. QR  4. Eigen  5. SVD  6. Diagonalization" << std::endl;
                    int sub; std::cin >> sub;
                    Matrix A = inputMatrix("Matrix A");
                    
                    if (sub == 1) {
                        auto [L, U] = Decomposer::LU(A);
                        std::cout << "L:" << std::endl; L.print();
                        std::cout << "U:" << std::endl; U.print();
                    } else if (sub == 2) {
                        Matrix L = Decomposer::Cholesky(A);
                        std::cout << "L:" << std::endl; L.print();
                    } else if (sub == 3) {
                        auto [Q, R] = Decomposer::QR(A);
                        std::cout << "Q:" << std::endl; Q.print();
                        std::cout << "R:" << std::endl; R.print();
                    } else if (sub == 4) {
                        std::cout << "1. Full Eigendecomposition (QR)  2. Dominant Eigenvalue (Power Iteration)" << std::endl;
                        int eigenSub; std::cin >> eigenSub;
                        if (eigenSub == 1) {
                            auto [Evals, Evecs] = Decomposer::Eigen(A);
                            std::cout << "Eigenvalues:" << std::endl; Evals.print();
                            std::cout << "Eigenvectors:" << std::endl; Evecs.print();
                        } else {
                            auto [lambda, v] = Decomposer::PowerIteration(A);
                            std::cout << "Dominant Eigenvalue: " << lambda << std::endl;
                            std::cout << "Dominant Eigenvector: "; VectorOps::print(v);
                        }
                    } else if (sub == 5) {
                        auto [U, S, V] = Decomposer::SVD(A);
                        std::cout << "U:" << std::endl; U.print();
                        std::cout << "S:" << std::endl; S.print();
                        std::cout << "V:" << std::endl; V.print();
                    } else if (sub == 6) {
                        auto [P, D] = Analysis::diagonalize(A);
                        std::cout << "P (Eigenvectors):" << std::endl; P.print();
                        std::cout << "D (Diagonal):" << std::endl; D.print();
                        std::cout << "Verify A = PDP^-1:" << std::endl;
                        try {
                            (P * D * LinearSolver::inverse(P)).print();
                        } catch (...) {
                            std::cout << "Could not verify (P might be singular)" << std::endl;
                        }
                    }
                    break;
                }
                case 6: { // Analysis
                    std::cout << "1. PCA  2. Rank  3. Trace  4. Condition Number  5. Quadratic Form" << std::endl;
                    int sub; std::cin >> sub;
                    if (sub == 1) {
                        Matrix X = inputMatrix("Data Matrix X (rows=samples, cols=features)");
                        int k;
                        std::cout << "Number of components: "; std::cin >> k;
                        auto [Comps, Transformed] = Analysis::PCA(X, k);
                        std::cout << "Principal Components:" << std::endl; Comps.print();
                        std::cout << "Transformed Data:" << std::endl; Transformed.print();
                    } else if (sub == 2) {
                        Matrix A = inputMatrix("Matrix A");
                        std::cout << "Rank: " << A.rank() << std::endl;
                    } else if (sub == 3) {
                        Matrix A = inputMatrix("Matrix A");
                        std::cout << "Trace: " << A.trace() << std::endl;
                    } else if (sub == 4) {
                        Matrix A = inputMatrix("Matrix A");
                        std::cout << "Condition Number: " << A.conditionNumber() << std::endl;
                    } else if (sub == 5) {
                        Matrix A = inputMatrix("Symmetric Matrix A");
                        std::vector<double> x = inputVector("Vector x");
                        try {
                            double val = QuadraticForm::evaluate(A, x);
                            std::cout << "Value (x^T A x): " << val << std::endl;
                            
                            auto def = QuadraticForm::analyzeDefiniteness(A);
                            std::cout << "Definiteness: " << QuadraticForm::definitenessToString(def) << std::endl;
                            
                            std::cout << "Canonical Form (Eigenvalues):" << std::endl;
                            QuadraticForm::canonicalForm(A).print();
                        } catch (const std::exception& e) {
                            std::cout << "Error: " << e.what() << std::endl;
                        }
                    }
                    break;
                }
                case 7: { // Statistics
                    std::cout << "1. Mean  2. Variance  3. Std Dev  4. Covariance Matrix  5. Correlation Matrix" << std::endl;
                    int sub; std::cin >> sub;
                    Matrix A = inputMatrix("Data Matrix (rows=samples, cols=features)");
                    
                    if (sub == 1) {
                        std::cout << "Mean (column-wise): ";
                        VectorOps::print(Statistics::mean(A, 0));
                    } else if (sub == 2) {
                        std::cout << "Variance (column-wise): ";
                        VectorOps::print(Statistics::variance(A, 0));
                    } else if (sub == 3) {
                        std::cout << "Standard Deviation (column-wise): ";
                        VectorOps::print(Statistics::stdDev(A, 0));
                    } else if (sub == 4) {
                        std::cout << "Covariance Matrix:" << std::endl;
                        Statistics::covarianceMatrix(A).print();
                    } else if (sub == 5) {
                        std::cout << "Correlation Matrix:" << std::endl;
                        Statistics::correlationMatrix(A).print();
                    }
                    break;
                }
                case 8: { // Matrix Manipulation
                    std::cout << "1. Submatrix  2. Horizontal Stack  3. Vertical Stack  4. Apply Function" << std::endl;
                    int sub; std::cin >> sub;
                    if (sub == 1) {
                        Matrix A = inputMatrix("Matrix A");
                        int sr, sc, nr, nc;
                        std::cout << "Start row, Start col, Num rows, Num cols: ";
                        std::cin >> sr >> sc >> nr >> nc;
                        Matrix Sub = A.submatrix(sr, sc, nr, nc);
                        std::cout << "Submatrix:" << std::endl;
                        Sub.print();
                    } else if (sub == 2) {
                        Matrix A = inputMatrix("Matrix A");
                        Matrix B = inputMatrix("Matrix B");
                        Matrix C = Matrix::hstack(A, B);
                        std::cout << "Horizontally stacked:" << std::endl;
                        C.print();
                    } else if (sub == 3) {
                        Matrix A = inputMatrix("Matrix A");
                        Matrix B = inputMatrix("Matrix B");
                        Matrix C = Matrix::vstack(A, B);
                        std::cout << "Vertically stacked:" << std::endl;
                        C.print();
                    } else if (sub == 4) {
                        Matrix A = inputMatrix("Matrix A");
                        std::cout << "Functions: 1. sin  2. cos  3. exp  4. square  5. sqrt  6. abs" << std::endl;
                        int func; std::cin >> func;
                        if (func == 1) {
                            Matrix Result = A.applyFunction([](double x) { return std::sin(x); });
                            std::cout << "Result:" << std::endl;
                            Result.print();
                        } else if (func == 2) {
                            Matrix Result = A.applyFunction([](double x) { return std::cos(x); });
                            std::cout << "Result:" << std::endl;
                            Result.print();
                        } else if (func == 3) {
                            Matrix Result = A.applyFunction([](double x) { return std::exp(x); });
                            std::cout << "Result:" << std::endl;
                            Result.print();
                        } else if (func == 4) {
                            Matrix Result = A.applyFunction([](double x) { return x * x; });
                            std::cout << "Result:" << std::endl;
                            Result.print();
                        } else if (func == 5) {
                            Matrix Result = A.applyFunction([](double x) { return std::sqrt(std::abs(x)); });
                            std::cout << "Result:" << std::endl;
                            Result.print();
                        } else if (func == 6) {
                            Matrix Result = A.applyFunction([](double x) { return std::abs(x); });
                            std::cout << "Result:" << std::endl;
                            Result.print();
                        }
                    }
                    break;
                }
                case 10: { // Sparse Matrices
                    std::cout << "1. Convert Dense to Sparse  2. Sparse Operations  3. View Statistics" << std::endl;
                    int sub; std::cin >> sub;
                    if (sub == 1) {
                        Matrix A = inputMatrix("Matrix A");
                        std::cout << "Convert to: 1. CSR  2. CSC" << std::endl;
                        int fmt; std::cin >> fmt;
                        if (fmt == 1) {
                            SparseMatrixCSR sparse = SparseMatrixCSR::fromDense(A);
                            sparse.printStats();
                            std::cout << "\nSparse matrix:" << std::endl;
                            sparse.print();
                        } else {
                            SparseMatrixCSC sparse = SparseMatrixCSC::fromDense(A);
                            sparse.printStats();
                            std::cout << "\nSparse matrix:" << std::endl;
                            sparse.print();
                        }
                    } else if (sub == 2) {
                        Matrix A = inputMatrix("Matrix A (will convert to sparse)");
                        SparseMatrixCSR sparseA = SparseMatrixCSR::fromDense(A);
                        std::cout << "Choose operation: 1. Scalar multiply  2. Add to another" << std::endl;
                        int op; std::cin >> op;
                        if (op == 1) {
                            double scalar;
                            std::cout << "Enter scalar: ";
                            std::cin >> scalar;
                            SparseMatrixCSR result = sparseA * scalar;
                            result.print();
                        } else {
                            Matrix B = inputMatrix("Matrix B");
                            SparseMatrixCSR sparseB = SparseMatrixCSR::fromDense(B);
                            SparseMatrixCSR result = sparseA + sparseB;
                            result.print();
                        }
                    } else if (sub == 3) {
                        Matrix A = inputMatrix("Matrix A");
                        SparseMatrixCSR sparse = SparseMatrixCSR::fromDense(A);
                        sparse.printStats();
                    }
                    break;
                }
                case 11: { // Iterative Solvers
                    std::cout << "1. Conjugate Gradient (CG)  2. GMRES  3. Preconditioned CG" << std::endl;
                    int sub; std::cin >> sub;
                    Matrix A = inputMatrix("Matrix A");
                    std::vector<double> b = inputVector("Vector b");
                    double tol;
                    int maxIter;
                    std::cout << "Enter tolerance (e.g., 1e-6): ";
                    std::cin >> tol;
                    std::cout << "Enter max iterations: ";
                    std::cin >> maxIter;
                    
                    std::vector<double> x;
                    if (sub == 1) {
                        std::cout << "Use sparse format? (1=Yes, 0=No): ";
                        int useSparse; std::cin >> useSparse;
                        if (useSparse) {
                            SparseMatrixCSR sparseA = SparseMatrixCSR::fromDense(A);
                            x = IterativeSolvers::CG(sparseA, b, tol, maxIter);
                        } else {
                            x = IterativeSolvers::CG(A, b, tol, maxIter);
                        }
                    } else if (sub == 2) {
                        int restart;
                        std::cout << "Enter restart parameter (e.g., 30): ";
                        std::cin >> restart;
                        std::cout << "Use sparse format? (1=Yes, 0=No): ";
                        int useSparse; std::cin >> useSparse;
                        if (useSparse) {
                            SparseMatrixCSR sparseA = SparseMatrixCSR::fromDense(A);
                            x = IterativeSolvers::GMRES(sparseA, b, restart, tol, maxIter);
                        } else {
                            x = IterativeSolvers::GMRES(A, b, restart, tol, maxIter);
                        }
                    } else if (sub == 3) {
                        x = IterativeSolvers::PCG(A, b, tol, maxIter);
                    }
                    
                    std::cout << "Solution x: ";
                    VectorOps::print(x);
                    break;
                }
                case 12: { // Complex Matrices
                    std::cout << "1. Create and Display  2. Arithmetic  3. Check Hermitian" << std::endl;
                    int sub; std::cin >> sub;
                    if (sub == 1) {
                        int r, c;
                        std::cout << "Enter dimensions (rows cols): ";
                        std::cin >> r >> c;
                        ComplexMatrix M(r, c);
                        std::cout << "Enter elements (real imag pairs):" << std::endl;
                        for (int i = 0; i < r; ++i) {
                            for (int j = 0; j < c; ++j) {
                                double re, im;
                                std::cin >> re >> im;
                                M(i, j) = std::complex<double>(re, im);
                            }
                        }
                        std::cout << "\nComplex Matrix:" << std::endl;
                        M.print();
                        std::cout << "\nConjugate Transpose:" << std::endl;
                        M.conjugateTranspose().print();
                    } else if (sub == 2) {
                        int r, c;
                        std::cout << "Enter dimensions for Matrix A (rows cols): ";
                        std::cin >> r >> c;
                        ComplexMatrix A(r, c);
                        std::cout << "Enter elements for A (real imag pairs):" << std::endl;
                        for (int i = 0; i < r; ++i) {
                            for (int j = 0; j < c; ++j) {
                                double re, im;
                                std::cin >> re >> im;
                                A(i, j) = std::complex<double>(re, im);
                            }
                        }
                        std::cout << "Operation: 1. Scalar multiply  2. Add to another" << std::endl;
                        int op; std::cin >> op;
                        if (op == 1) {
                            double re, im;
                            std::cout << "Enter scalar (real imag): ";
                            std::cin >> re >> im;
                            ComplexMatrix result = A * std::complex<double>(re, im);
                            result.print();
                        } else {
                            int r2, c2;
                            std::cout << "Enter dimensions for Matrix B (rows cols): ";
                            std::cin >> r2 >> c2;
                            ComplexMatrix B(r2, c2);
                            std::cout << "Enter elements for B (real imag pairs):" << std::endl;
                            for (int i = 0; i < r2; ++i) {
                                for (int j = 0; j < c2; ++j) {
                                    double re, im;
                                    std::cin >> re >> im;
                                    B(i, j) = std::complex<double>(re, im);
                                }
                            }
                            ComplexMatrix result = A + B;
                            result.print();
                        }
                    } else if (sub == 3) {
                        int n;
                        std::cout << "Enter size for square matrix: ";
                        std::cin >> n;
                        ComplexMatrix M(n, n);
                        std::cout << "Enter elements (real imag pairs):" << std::endl;
                        for (int i = 0; i < n; ++i) {
                            for (int j = 0; j < n; ++j) {
                                double re, im;
                                std::cin >> re >> im;
                                M(i, j) = std::complex<double>(re, im);
                            }
                        }
                        M.print();
                        if (M.isHermitian()) {
                            std::cout << "Matrix is Hermitian!" << std::endl;
                        } else {
                            std::cout << "Matrix is NOT Hermitian." << std::endl;
                        }
                    }
                    break;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }
    return 0;
}
