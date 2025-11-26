# Linear Algebra CLI Tool

A comprehensive command-line tool for linear algebra operations with step-by-step explanations, supporting matrix arithmetic, decompositions, eigenvalue computations, and statistical analysis.

## Features

### Matrix Operations
- **Basic Arithmetic**: Addition, Subtraction, Multiplication, Scalar Multiplication
- **Transpose**: Matrix transposition
- **Power**: Fast exponentiation using binary exponentiation
- **Rank & Trace**: Compute matrix rank and trace

### Vector Operations
- **Dot Product**: Inner product of two vectors
- **Norm**: Euclidean norm (L2 norm)
- **Cross Product**: Cross product for 3D vectors
- **Normalization**: Unit vector computation

### Linear Systems
- **Gaussian Elimination**: Solve $Ax = b$ with partial pivoting
- **Determinant**: Compute determinant using Gaussian elimination
- **Matrix Inverse**: Find inverse using Gauss-Jordan elimination

### Decompositions
- **LU Decomposition**: $A = LU$ factorization
- **Cholesky Decomposition**: $A = LL^T$ for positive-definite matrices
- **QR Decomposition**: $A = QR$ using Gram-Schmidt orthogonalization
- **Eigendecomposition**: 
  - Full eigendecomposition using QR algorithm (handles both symmetric and non-symmetric matrices)
  - Power iteration for dominant eigenvalue/eigenvector
- **SVD**: Singular Value Decomposition $A = U\Sigma V^T$
- **Diagonalization**: Express $A = PDP^{-1}$ where $D$ is diagonal

### Analysis & Statistics
- **PCA**: Principal Component Analysis for dimensionality reduction
- **Rank**: Matrix rank computation
- **Trace**: Sum of diagonal elements

### Verbose Mode
- **Step-by-Step Explanations**: Toggle verbose mode (Option 7) to see detailed algorithmic steps
- Logs row operations during Gaussian elimination
- Shows convergence progress for iterative methods
- Displays intermediate results for decompositions

## Building

To build the project, simply run:

```bash
make
```

This compiles all source files and creates the `linear_algebra` executable.

To clean build artifacts:

```bash
make clean
```

## Running

Start the program with:

```bash
./linear_algebra
```

## Usage Examples

### Example 1: Solving a Linear System

```
Select option: 3
Enter dimensions for Matrix A (rows cols): 2 2
Enter elements row by row:
2 1
1 3
Enter size for Vector b: 2
Enter elements: 4 7
Solution x: ( 1 2 )
```

### Example 2: Eigendecomposition with Verbose Mode

```
Select option: 7  # Toggle verbose mode ON
Select option: 5  # Decompositions
1. LU  2. Cholesky  3. QR  4. Eigen  5. SVD  6. Diagonalization
4
Enter dimensions for Matrix A (rows cols): 2 2
Enter elements row by row:
4 1
2 3

[STEP] Starting Eigen Decomposition (QR Algorithm)
[STEP] Starting QR Decomposition (Gram-Schmidt)
...
[STEP] Converged at iteration 5
[STEP] Matrix is symmetric, Q contains eigenvectors
Eigenvalues:
[ 5.0000 0.0000 ]
[ 0.0000 2.0000 ]
Eigenvectors:
[ 0.8944 -0.4472 ]
[ 0.4472 0.8944 ]
```

### Example 3: Power Iteration

```
Select option: 5
4
1. Full Eigendecomposition (QR)  2. Dominant Eigenvalue (Power Iteration)
2
Enter dimensions for Matrix A (rows cols): 2 2
Enter elements row by row:
4 1
2 3
Dominant Eigenvalue: 5
Dominant Eigenvector: ( 0.8944 0.4472 )
```

## Technical Details

### Algorithms
- **Gaussian Elimination**: Partial pivoting for numerical stability
- **QR Algorithm**: Iterative Schur decomposition with back-substitution for eigenvectors
- **Power Iteration**: Rayleigh quotient for dominant eigenvalue estimation
- **Gram-Schmidt**: Classical orthogonalization for QR decomposition

### Numerical Precision
- Tolerance: `1e-10` for convergence checks
- Pivot threshold: `1e-10` for singularity detection
- Maximum iterations: `1000` for iterative methods

### Matrix Support
- **Symmetric matrices**: Optimized eigenvector computation
- **Non-symmetric matrices**: Back-substitution from Schur form
- **Positive-definite matrices**: Cholesky decomposition
- **Rectangular matrices**: SVD and QR decomposition

## Project Structure

```
.
├── main.cpp           # CLI interface and menu system
├── Matrix.h/cpp       # Matrix class with basic operations
├── VectorOps.h        # Vector operations (header-only)
├── LinearSolver.h/cpp # Linear systems and matrix inverse
├── Decomposer.h/cpp   # Matrix decompositions
├── Analysis.h/cpp     # PCA, diagonalization, transformations
├── Logger.h           # Step-by-step logging (header-only)
├── Makefile           # Build configuration
└── README.md          # This file
```

## Notes

- For best results with eigendecomposition, symmetric matrices are recommended
- Power iteration finds the eigenvalue with largest absolute value
- Cholesky decomposition requires symmetric positive-definite matrices
- Use verbose mode (Option 7) to understand algorithm behavior

