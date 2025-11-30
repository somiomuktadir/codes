#include "../include/Matrix.h"
#include "TestFramework.h"

using namespace LinAlg;

TEST(Matrix_Construction) {
    Matrix A(3, 3, 1.0);
    ASSERT_EQ(A.getRows(), 3);
    ASSERT_EQ(A.getCols(), 3);
    ASSERT_EQ(A(0, 0), 1.0);
    ASSERT_EQ(A(2, 2), 1.0);

    Matrix B = {{1, 2}, {3, 4}};
    ASSERT_EQ(B.getRows(), 2);
    ASSERT_EQ(B.getCols(), 2);
    ASSERT_EQ(B(0, 1), 2.0);
    ASSERT_EQ(B(1, 0), 3.0);
}

TEST(Matrix_Arithmetic) {
    Matrix A = {{1, 2}, {3, 4}};
    Matrix B = {{5, 6}, {7, 8}};
    
    Matrix C = A + B;
    ASSERT_EQ(C(0, 0), 6.0);
    ASSERT_EQ(C(1, 1), 12.0);
    
    Matrix D = A * B;
    // [1 2] * [5 6] = [19 22]
    // [3 4]   [7 8]   [43 50]
    ASSERT_EQ(D(0, 0), 19.0);
    ASSERT_EQ(D(0, 1), 22.0);
    ASSERT_EQ(D(1, 0), 43.0);
    ASSERT_EQ(D(1, 1), 50.0);
}

TEST(Matrix_Transpose) {
    Matrix A = {{1, 2, 3}, {4, 5, 6}};
    Matrix T = A.transpose();
    
    ASSERT_EQ(T.getRows(), 3);
    ASSERT_EQ(T.getCols(), 2);
    ASSERT_EQ(T(0, 1), 4.0);
    ASSERT_EQ(T(2, 0), 3.0);
}

TEST(Matrix_Norms) {
    Matrix A = {{1, -2}, {-3, 4}};
    // Frobenius: sqrt(1+4+9+16) = sqrt(30)
    ASSERT_NEAR(A.frobeniusNorm(), std::sqrt(30.0), 1e-10);
}
