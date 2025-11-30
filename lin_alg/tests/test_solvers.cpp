#include "../include/LinearSolver.h"
#include "../include/Matrix.h"
#include "TestFramework.h"

using namespace LinAlg;

TEST(LinearSolver_PseudoInverse) {
    // A = [1 0; 0 1; 0 0]
    Matrix A(3, 2);
    A(0, 0) = 1.0; A(1, 1) = 1.0;
    
    // A^+ should be [1 0 0; 0 1 0]
    Matrix P = LinearSolver::pseudoInverse(A);
    
    ASSERT_EQ(P.getRows(), 2);
    ASSERT_EQ(P.getCols(), 3);
    ASSERT_NEAR(P(0, 0), 1.0, 1e-9);
    ASSERT_NEAR(P(1, 1), 1.0, 1e-9);
    ASSERT_NEAR(P(0, 2), 0.0, 1e-9);
}

TEST(LinearSolver_SolvePseudoInverse) {
    // Overdetermined system
    // x = 1, y = 1
    // x + 0y = 1
    // 0x + y = 1
    // x + y = 2
    
    Matrix A(3, 2);
    A(0, 0) = 1.0; A(0, 1) = 0.0;
    A(1, 0) = 0.0; A(1, 1) = 1.0;
    A(2, 0) = 1.0; A(2, 1) = 1.0;
    
    std::vector<double> b = {1.0, 1.0, 2.0};
    
    std::vector<double> x = LinearSolver::solvePseudoInverse(A, b);
    
    ASSERT_NEAR(x[0], 1.0, 1e-9);
    ASSERT_NEAR(x[1], 1.0, 1e-9);
}
