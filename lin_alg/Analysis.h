#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "Matrix.h"
#include <vector>
#include <utility>

class Analysis {
public:
    // Principal Component Analysis
    // Input: Data matrix X (rows = samples, cols = features)
    // Returns: {Principal Components (Eigenvectors), Transformed Data}
    static std::pair<Matrix, Matrix> PCA(const Matrix& X, int numComponents);
    
    // Returns {P, D} such that A = P * D * P^-1
    static std::pair<Matrix, Matrix> diagonalize(const Matrix& A);

    // Apply Linear Transformation
    // Input: Matrix A (transformation), Vector v
    // Returns: Transformed vector
    static std::vector<double> transform(const Matrix& A, const std::vector<double>& v);
    
    // Apply Linear Transformation to a set of points (Matrix)
    static Matrix transform(const Matrix& T, const Matrix& Points);
};

#endif // ANALYSIS_H
