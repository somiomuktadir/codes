#include "QuadraticForm.h"
#include "Decomposer.h"
#include "VectorOps.h"
#include <stdexcept>
#include <algorithm>

namespace LinAlg {

double QuadraticForm::evaluate(const Matrix& A, const std::vector<double>& x) {
    if (!A.isSquare()) throw std::invalid_argument("Matrix must be square for quadratic form");
    if (A.getRows() != (int)x.size()) throw std::invalid_argument("Dimension mismatch between matrix and vector");

    // x^T * A * x
    // Let y = A * x
    std::vector<double> y(x.size(), 0.0);
    for (int i = 0; i < A.getRows(); ++i) {
        for (int j = 0; j < A.getCols(); ++j) {
            y[i] += A(i, j) * x[j];
        }
    }
    
    return VectorOps::dot(x, y);
}

QuadraticForm::Definiteness QuadraticForm::analyzeDefiniteness(const Matrix& A) {
    if (!A.isSymmetric()) throw std::invalid_argument("Matrix must be symmetric for definiteness analysis");
    
    // Use Eigenvalues
    auto [Evals, Evecs] = Decomposer::Eigen(A);
    
    bool allPos = true;
    bool allNeg = true;
    bool hasZero = false;
    
    for (int i = 0; i < Evals.getRows(); ++i) {
        double lambda = Evals(i, i);
        if (std::abs(lambda) < Matrix::epsilon()) {
            hasZero = true;
        } else if (lambda > 0) {
            allNeg = false;
        } else {
            allPos = false;
        }
    }
    
    if (allPos && !hasZero) return Definiteness::POSITIVE_DEFINITE;
    if (allPos && hasZero) return Definiteness::POSITIVE_SEMIDEFINITE;
    if (allNeg && !hasZero) return Definiteness::NEGATIVE_DEFINITE;
    if (allNeg && hasZero) return Definiteness::NEGATIVE_SEMIDEFINITE;
    
    return Definiteness::INDEFINITE;
}

std::string QuadraticForm::definitenessToString(Definiteness d) {
    switch (d) {
        case Definiteness::POSITIVE_DEFINITE: return "Positive Definite";
        case Definiteness::POSITIVE_SEMIDEFINITE: return "Positive Semi-Definite";
        case Definiteness::NEGATIVE_DEFINITE: return "Negative Definite";
        case Definiteness::NEGATIVE_SEMIDEFINITE: return "Negative Semi-Definite";
        case Definiteness::INDEFINITE: return "Indefinite";
        default: return "Unknown";
    }
}

Matrix QuadraticForm::canonicalForm(const Matrix& A) {
    if (!A.isSymmetric()) throw std::invalid_argument("Matrix must be symmetric for canonical form");
    
    // The canonical form is the diagonal matrix of eigenvalues
    auto [Evals, Evecs] = Decomposer::Eigen(A);
    return Evals;
}

} // namespace LinAlg
