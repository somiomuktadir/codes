#include "../include/QuadraticForm.h"
#include "../include/Matrix.h"
#include "TestFramework.h"

using namespace LinAlg;

TEST(QuadraticForm_Evaluate) {
    // A = I
    Matrix A = Matrix::identity(2);
    std::vector<double> x = {1.0, 2.0};
    
    // x^T I x = 1*1 + 2*2 = 5
    double val = QuadraticForm::evaluate(A, x);
    ASSERT_NEAR(val, 5.0, 1e-9);
}

TEST(QuadraticForm_Definiteness) {
    // Positive Definite: [2 0; 0 2]
    Matrix Pos(2, 2);
    Pos(0, 0) = 2.0; Pos(1, 1) = 2.0;
    ASSERT_TRUE(QuadraticForm::analyzeDefiniteness(Pos) == QuadraticForm::Definiteness::POSITIVE_DEFINITE);
    
    // Negative Definite: [-2 0; 0 -2]
    Matrix Neg(2, 2);
    Neg(0, 0) = -2.0; Neg(1, 1) = -2.0;
    ASSERT_TRUE(QuadraticForm::analyzeDefiniteness(Neg) == QuadraticForm::Definiteness::NEGATIVE_DEFINITE);
    
    // Indefinite: [2 0; 0 -2]
    Matrix Ind(2, 2);
    Ind(0, 0) = 2.0; Ind(1, 1) = -2.0;
    ASSERT_TRUE(QuadraticForm::analyzeDefiniteness(Ind) == QuadraticForm::Definiteness::INDEFINITE);
}
