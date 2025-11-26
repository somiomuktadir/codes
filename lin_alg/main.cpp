#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include "Matrix.h"
#include "VectorOps.h"
#include "LinearSolver.h"
#include "Decomposer.h"
#include "Decomposer.h"
#include "Analysis.h"
#include "Logger.h"

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
            std::cin >> M.at(i, j);
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
    std::cout << "\n=== Linear Algebra CLI ===" << std::endl;
    std::cout << "1. Matrix Operations (+, -, *, T)" << std::endl;
    std::cout << "2. Vector Operations (Dot, Norm, Cross)" << std::endl;
    std::cout << "3. Linear Solver (Ax = b)" << std::endl;
    std::cout << "4. Determinant & Inverse" << std::endl;
    std::cout << "5. Decompositions (LU, Cholesky, QR, Eigen, SVD, Diagonalization)" << std::endl;
    std::cout << "6. Analysis (PCA, Rank, Trace)" << std::endl;
    std::cout << "7. Toggle Verbose Mode (Current: " << (Logger::getInstance().isEnabled() ? "ON" : "OFF") << ")" << std::endl;
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
            
            if (choice == 7) {
                if (Logger::getInstance().isEnabled()) Logger::getInstance().disable();
                else Logger::getInstance().enable();
                continue;
            }
            
            switch (choice) {
                case 1: { // Matrix Ops
                    std::cout << "1. Add  2. Subtract  3. Multiply  4. Transpose  5. Power" << std::endl;
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
                    Matrix A = inputMatrix("Matrix A");
                    std::vector<double> b = inputVector("Vector b");
                    std::vector<double> x = LinearSolver::solve(A, b);
                    std::cout << "Solution x: ";
                    VectorOps::print(x);
                    break;
                }
                case 4: { // Det & Inv
                    Matrix A = inputMatrix("Matrix A");
                    std::cout << "Determinant: " << LinearSolver::determinant(A) << std::endl;
                    try {
                        std::cout << "Inverse:" << std::endl;
                        LinearSolver::inverse(A).print();
                    } catch (const std::exception& e) {
                        std::cout << "Inverse not possible: " << e.what() << std::endl;
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
                    std::cout << "1. PCA  2. Rank  3. Trace" << std::endl;
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
