# Linear Algebra CLI Tool

A high-performance, comprehensive C++ library and command-line interface for advanced linear algebra operations. This tool is designed for educational and practical purposes, offering step-by-step algorithmic explanations alongside robust numerical computations.

## Overview

The Linear Algebra CLI Tool provides a suite of operations ranging from basic matrix arithmetic to complex decompositions and iterative solvers. It is built with performance in mind, utilizing contiguous memory storage, cache-friendly access patterns, and optional BLAS/OpenMP integration.

## Key Features

### Core Matrix Operations
-   **Arithmetic**: Addition, subtraction, multiplication, and scalar operations.
-   **Element-wise Operations**: Hadamard product and function application (sin, cos, exp, etc.).
-   **Manipulation**: Submatrix extraction, horizontal/vertical stacking, and resizing.
-   **Vector Operations**: Dot product, L2 norm, cross product (3D), and projection.

### Advanced Decompositions
-   **LU Decomposition**: $PA = LU$ factorization with partial pivoting.
-   **Cholesky Decomposition**: $A = LL^T$ for symmetric positive-definite matrices.
-   **QR Decomposition**: $A = QR$ using Modified Gram-Schmidt orthogonalization.
-   **Eigendecomposition**:
    -   QR Algorithm with Hessenberg reduction for full spectrum analysis.
    -   Power Iteration for dominant eigenvalue estimation.
    -   Schur Decomposition for upper triangular forms.
-   **Singular Value Decomposition (SVD)**: $A = U\Sigma V^T$ via Golub-Kahan Bidiagonalization.
-   **Diagonalization**: $A = PDP^{-1}$ spectral decomposition.

### Linear System Solvers
-   **Direct Solvers**:
    -   Gaussian Elimination with partial pivoting.
    -   Cholesky Solver for symmetric positive-definite systems.
    -   Least Squares solver for overdetermined systems via QR.
-   **Iterative Solvers**:
    -   **Conjugate Gradient (CG)**: Optimized for symmetric positive-definite systems.
    -   **GMRES**: Generalized Minimal Residual method for non-symmetric systems.
    -   **Preconditioning**: Jacobi preconditioner support for improved convergence.

### Sparse Matrix Support
Efficient handling of large, sparse datasets using compressed storage formats:
-   **CSR (Compressed Sparse Row)**: Optimized for row-wise operations and matrix-vector multiplication.
-   **CSC (Compressed Sparse Column)**: Optimized for column-wise operations.
-   **Operations**: Memory-efficient arithmetic and seamless conversion between dense and sparse formats.
-   **I/O**: Support for the NIST Matrix Market (`.mtx`) format.

### Complex Number Support
Full support for complex-valued matrices:
-   **Arithmetic**: Complex addition, multiplication, and scalar operations.
-   **Hermitian Operations**: Conjugate transpose and Hermitian property verification.

### Data Analysis & Statistics
-   **Principal Component Analysis (PCA)**: Dimensionality reduction and feature extraction.
-   **Matrix Statistics**: Computation of rank, trace, and condition number.

## Technical Highlights

-   **Performance**:
    -   **Flat Storage**: Uses `std::vector<double>` for contiguous memory layout, improving cache locality.
    -   **Parallelization**: Optional OpenMP support for multi-threaded operations on large matrices.
    -   **BLAS Integration**: Optional linkage with BLAS libraries for accelerated matrix multiplication.
-   **Numerical Stability**:
    -   Uses Modified Gram-Schmidt for QR decomposition.
    -   Implements Golub-Kahan Bidiagonalization for SVD.
    -   Includes iterative refinement for linear system solutions.
-   **Verbose Mode**: A unique educational feature that logs detailed, step-by-step algorithmic progress (e.g., row operations, convergence steps) to the console.

## Building the Project

The project uses a standard `Makefile` for compilation.

### Prerequisites
-   A C++17 compatible compiler (GCC, Clang).
-   (Optional) OpenMP and BLAS for performance enhancements.

### Compilation

To build the release version with optimizations:

```bash
make full
```

To clean build artifacts:

```bash
make clean
```

## Usage

Run the interactive command-line interface:

```bash
./linear_algebra
```

### CLI Menu Structure

The interactive menu provides access to all features:

1.  **Matrix Operations**: Basic arithmetic, transpose, power.
2.  **Vector Operations**: Dot/cross products, norms.
3.  **Linear Solver**: Direct solvers (Gaussian, Cholesky, Least Squares).
4.  **Determinant & Inverse**: Matrix inversion and determinant calculation.
5.  **Decompositions**: LU, QR, Cholesky, Eigen, SVD.
6.  **Analysis**: PCA, Rank, Condition Number.
7.  **File I/O**: CSV and Matrix Market import/export.
8.  **Matrix Manipulation**: Resizing, stacking, submatrices.
9.  **Toggle Verbose Mode**: Enable/disable educational logging.
10. **Sparse Matrices**: CSR/CSC operations and statistics.
11. **Iterative Solvers**: CG and GMRES methods.
12. **Complex Matrices**: Complex arithmetic and Hermitian checks.

### Example: Solving a Linear System

```bash
Select option: 3
1. Standard  2. Cholesky (for SPD)  3. Least Squares  4. Refined
1
Enter dimensions for Matrix A (rows cols): 2 2
Enter elements row by row:
2 1
1 3
Enter size for Vector b: 2
Enter elements: 4 7
Solution x: ( 1.0000 2.0000 )
```

## API Reference

The library is encapsulated within the `LinAlg` namespace.

```cpp
#include "Matrix.h"
#include "LinearSolver.h"

using namespace LinAlg;

int main() {
    // Create a 3x3 matrix
    Matrix A = {{4, 1, 1}, {1, 4, 1}, {1, 1, 4}};
    std::vector<double> b = {6, 6, 6};

    // Solve Ax = b
    std::vector<double> x = LinearSolver::solve(A, b);

    // Print result
    VectorOps::print(x);
    return 0;
}
```

## File Formats

-   **CSV**: Standard comma-separated values for dense matrices.
-   **Matrix Market (.mtx)**: Industry-standard format for sparse and dense matrices, supporting coordinate and array formats.