#ifndef QUADRATIC_FORM_H
#define QUADRATIC_FORM_H

#include "Matrix.h"
#include <vector>
#include <string>

namespace LinAlg {

class QuadraticForm {
public:
    // Compute x^T * A * x
    static double evaluate(const Matrix& A, const std::vector<double>& x);

    enum class Definiteness {
        POSITIVE_DEFINITE,
        POSITIVE_SEMIDEFINITE,
        NEGATIVE_DEFINITE,
        NEGATIVE_SEMIDEFINITE,
        INDEFINITE,
        UNKNOWN
    };

    // Determine definiteness using eigenvalues
    static Definiteness analyzeDefiniteness(const Matrix& A);
    
    static std::string definitenessToString(Definiteness d);

    // Returns the diagonal matrix D in the canonical form Q = y^T D y
    // where y = P^T x, and P is the orthogonal matrix of eigenvectors
    static Matrix canonicalForm(const Matrix& A);
};

} // namespace LinAlg

#endif // QUADRATIC_FORM_H
