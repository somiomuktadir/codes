#include "Analysis.h"
#include "Decomposer.h"
#include "VectorOps.h"
#include <numeric>
#include <algorithm>

std::pair<Matrix, Matrix> Analysis::PCA(const Matrix& X, int numComponents) {
    int n = X.getRows();
    int m = X.getCols();
    
    // 1. Center the data (subtract mean from each column)
    Matrix Centered = X;
    for (int j = 0; j < m; ++j) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += X.at(i, j);
        }
        double mean = sum / n;
        for (int i = 0; i < n; ++i) {
            Centered.at(i, j) -= mean;
        }
    }
    
    // 2. Compute Covariance Matrix: C = (X^T * X) / (n - 1)
    Matrix Cov = (Centered.transpose() * Centered) * (1.0 / (n - 1));
    
    // 3. Eigen Decomposition of Covariance Matrix
    auto [Evals, Evecs] = Decomposer::Eigen(Cov);
    
    // 4. Sort Eigenvectors by Eigenvalues (descending)
    // Note: Decomposer::Eigen returns unsorted. We need to sort.
    std::vector<std::pair<double, int>> evalPairs;
    for (int i = 0; i < m; ++i) {
        evalPairs.push_back({Evals.at(i, i), i});
    }
    std::sort(evalPairs.rbegin(), evalPairs.rend());
    
    Matrix SortedEvecs(m, m);
    for (int i = 0; i < m; ++i) {
        int oldIndex = evalPairs[i].second;
        for (int k = 0; k < m; ++k) {
            SortedEvecs.at(k, i) = Evecs.at(k, oldIndex);
        }
    }
    
    // 5. Select top 'numComponents' eigenvectors
    Matrix Components(m, numComponents);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < numComponents; ++j) {
            Components.at(i, j) = SortedEvecs.at(i, j);
        }
    }
    
    // 6. Project data: Y = X_centered * Components
    Matrix Transformed = Centered * Components;
    
    return {Components, Transformed};
}

std::vector<double> Analysis::transform(const Matrix& A, const std::vector<double>& v) {
    int rows = A.getRows();
    int cols = A.getCols();
    if (v.size() != cols) throw std::invalid_argument("Vector dimension mismatch for transformation");
    
    std::vector<double> result(rows, 0.0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i] += A.at(i, j) * v[j];
        }
    }
    return result;
}

Matrix Analysis::transform(const Matrix& T, const Matrix& Points) {
    // Points are assumed to be row vectors in the Points matrix? 
    // Usually T * v. If Points has points as columns, then T * Points.
    // If Points has points as rows, then Points * T^T.
    // Let's assume Points are columns for standard linear algebra notation T * P.
    return T * Points;
}

std::pair<Matrix, Matrix> Analysis::diagonalize(const Matrix& A) {
    // A = P * D * P^-1
    // Decomposer::Eigen returns {D, P} (Eigenvalues, Eigenvectors)
    // We want to return {P, D}
    auto [D, P] = Decomposer::Eigen(A);
    return {P, D};
}
